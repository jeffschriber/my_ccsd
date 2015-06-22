/*
 *@BEGIN LICENSE
 *
 * my_ccsd by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include <libdpd/dpd.h>
#include "psifiles.h"
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>
#include <libmints/petitelist.h>
#include "boost/multi_array.hpp"


// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)


INIT_PLUGIN

namespace psi{ namespace my_ccsd{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "MY_CCSD" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);

        //convergence criteria
        options.add_double("e_conv", 1e-10);
        options.add_double("t2_conv", 1e-6);

        //Lets user run calculation with (T)
        options.add_bool("P_TRIPLES", false);

        /*options for SRG-CCSD*/
        options.add_bool("SRG",false);
        options.add_double("srg_s", 1.0);
    }

    return true;
}


extern "C" PsiReturnType
my_ccsd(Options &options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");
   
    // Grab the global (default) PSIO object, for file I/O
    boost::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Now we want the reference (SCF) wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");

    // Quickly check that there are no open shell orbitals here...
    int nirrep  = wfn->nirrep();

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example in the test suite.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());



    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");

    //Make 1D array to store mo_TEIs
    int nmo = wfn->nmo();
    int dim = wfn->nmo();
    int dim2 = dim*dim*dim*dim;
    SharedVector TEI(new Vector("TEI array", dim2));

    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // Print out the absolute orbital numbers, the relative (within irrep)
                // numbers, the symmetries, and the integral itself
                //psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f, "
                                 "symmetries = (%1d %1d | %1d %1d), "
                                 "relative indices = (%2d %2d | %2d %2d)\n",
                                // p, q, r, s, K.matrix[h][pq][rs],
                                 //psym, qsym, rsym, ssym,
                                 //prel, qrel, rrel, srel);
                TEI->set(INDEX4(p,q,r,s),K.matrix[h][pq][rs]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    //Print some stuff
    outfile->Printf("\n     ======================================");
    outfile->Printf("\n     =                                    =");
    outfile->Printf("\n     =           A CCSD Plugin            =");
    outfile->Printf("\n     =      Written by: Jeff Schriber     =");
    outfile->Printf("\n     =                                    =");
    outfile->Printf("\n     =            Version 0.1             =");
    outfile->Printf("\n     ======================================\n\n");



    //Convert to spin orbital basis
    int pr;
    int qs;
    int prqs;
    int ps;
    int qr;
    int psqr;

    double value1;
    double value2;

    outfile->Flush();
    int nso = 2*dim;

    //Here I define 'quad', which creates a 4D array of doubles
    typedef boost::multi_array<double, 4> quad;
    quad spo_TEI(boost::extents[nso][nso][nso][nso]);

    //Construct the antisymmetrized TEIs in spin orbital basis
    for(int p=0; p < nso; p++)
      for(int q=0; q < nso; q++)
        for(int r=0; r < nso; r++)
          for(int s=0; s < nso; s++) {
            pr = INDEX2(p/2,r/2);
            qs = INDEX2(q/2,s/2);
            prqs = INDEX2(pr,qs);
            value1 = TEI->get(prqs) * (p%2 == r%2) * (q%2 == s%2);
            ps = INDEX2(p/2,s/2);
            qr = INDEX2(q/2,r/2);
            psqr = INDEX2(ps,qr);
            value2 = TEI->get(psqr) * (p%2 == s%2) * (q%2 == r%2);
            spo_TEI[p][q][r][s] = value1 - value2;
          }

    /*Now for the one-electron integrals, we need to convert
     * the Fock matrix to the spin orbital basis*/

    boost::shared_ptr<Matrix> cMat(new Matrix("CMat",nmo, nmo));
    boost::shared_ptr<Matrix> fMatAO(new Matrix("fMatAO", nmo, nmo));
    boost::shared_ptr<Matrix> fMat(new Matrix("fMat", 2*nmo, 2*nmo));
    boost::shared_ptr<Matrix> fMatMO(new Matrix("fMatMO", nmo, nmo));

    cMat = wfn->Ca();
    fMatAO = wfn->Fa();
    fMatMO->transform(cMat, fMatAO, cMat);

    //Map the Fock matrix from MO basis to spin orbital basis
    for(int i=0; i<nmo; i++){
        fMat->set(2*i,2*i,fMatMO->get(i,i));
        fMat->set(2*i+1,2*i+1,fMatMO->get(i,i));
    }

    //Generate the indices, nocc and nvir
    boost::shared_ptr<Vector> eps_occ = wfn->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wfn->epsilon_a_subset("AO", "VIR");

    int nocc = 2 * eps_occ->dim();
    int nvir = 2 * eps_vir->dim();


    //Build the initial guess amplitudes.
    boost::shared_ptr<Matrix> t1(new Matrix("t1", nocc, nvir));
    quad t2(boost::extents[nocc][nocc][nvir][nvir]);

    t1->zero(); //Initial guess for t1 amplitudes is zero


    //Generate the initial guess for t2 amplitudes
    for(int i=0; i<nocc; i++){
        for(int j=0; j<nocc; j++){
            for(int a = 0; a<nvir;a++ ){
                for(int b=0; b<nvir; b++){
                    t2[i][j][a][b] = (spo_TEI[i][j][a+nocc][b+nocc])/(fMat->get(i,i) + fMat->get(j,j) - fMat->get(a+nocc,a+nocc) - fMat->get(b+nocc, b+nocc));
                }
            }
        }
    }

    //Calculate the MP2 energy.
    double Emp2 = 0.0;
    for(int i=0; i<nocc; i++){
        for(int j=0; j<nocc; j++){
            for(int a = 0; a<nvir;a++ ){
                for(int b = 0; b<nvir; b++){
                   Emp2 += t2[i][j][a][b]*spo_TEI[i][j][a+nocc][b+nocc];
                }
            }
        }
    }

    Emp2 = Emp2/(4.0);

    //Grab RHF reference energy
    double Eref = wfn->reference_energy();


    bool converged = false;
    double Ecc = 0.0;
    double Eprev = 0.0;
    double dE = 0.0;
    int count = 0;
    quad t2prev(boost::extents[nocc][nocc][nvir][nvir]);

    if(options.get_bool("SRG") == true){
        outfile->Printf("\n Computing SRG Transformation\n");
    }

    outfile->Printf("\n    Iter                E(CCSD)                Delta(E)            RMS(T2)   \n");
    outfile->Printf("----------------------------------------------------------------------------------");
    outfile->Flush();

    //Begin the loop to generate updated t1 and t2 amps, break when converged
    while(converged == false){

        /*Definitions of the one-particle intermediates. Names are based on indices used in the Stanton paper.
         *All intermediates are build term-by-term, usually using several sets of loops.
        */

        boost::shared_ptr<Matrix> Fae(new Matrix("Fae", nvir, nvir));
        boost::shared_ptr<Matrix> Fmi(new Matrix("Fmi", nocc, nocc));
        boost::shared_ptr<Matrix> Fme(new Matrix("Fme", nocc, nvir));

        quad tau(boost::extents[nocc][nocc][nvir][nvir]);
        quad taut(boost::extents[nocc][nocc][nvir][nvir]);

        taut.empty();
        tau.empty();

        for(int i = 0; i < nocc; i++){
            for(int j = 0; j<  nocc; j++){
                for(int a = 0; a < nvir; a++ ){
                    for(int b = 0; b < nvir; b++){
                        taut[i][j][a][b] = t2[i][j][a][b] + (0.5)*t1->get(i,a)*t1->get(j,b) - (0.5)*t1->get(i,b)*t1->get(j,a);
                        tau[i][j][a][b] = t2[i][j][a][b] + t1->get(i,a)*t1->get(j,b) - t1->get(i,b)*t1->get(j,a);
                    }
                }
            }
        }

        //Calculate Fae

        int del = 0;
        for(int a = 0; a < nvir; a++){
            for(int e = 0; e < nvir; e++){
                if(a == e){ del = 1; }
                Fae->set(a,e, fMat->get(a+nocc,e+nocc)*(1-del));
            }
        }

        for(int m = 0; m < nocc; m++){
            for(int a = 0; a < nvir; a++){
                for(int e = 0; e < nvir; e++){
                    Fae->add(a, e, (-0.5)*fMat->get(m,e+nocc)*t1->get(m,a));
                }
            }
        }

        for(int m = 0; m < nocc; m++){
            for(int f = 0; f < nvir; f++){
                for(int a = 0; a < nvir; a++){
                    for(int e = 0; e < nvir; e++){
                        Fae->add(a,e, t1->get(m,f)*spo_TEI[m][a+nocc][f+nocc][e+nocc]);
                    }
                }
             }
        }

        for(int m = 0; m < nocc; m++){
            for(int n = 0; n < nocc; n++){
                for(int f = 0; f < nvir; f++){
                    for(int a = 0; a < nvir; a++){
                        for(int e = 0; e < nvir; e++){
                            Fae->add(a,e, (-0.5)*taut[m][n][a][f]*spo_TEI[m][n][e+nocc][f+nocc]);
                        }
                    }
                }
            }
        }

        //Build Fmi
        del = 0;
        for(int m = 0; m < nocc; m++){
            for(int i = 0; i < nocc; i++){
                if(m == i){ del = 1; }
                Fmi->set(m,i, fMat->get(m,i)*(1-del));
            }
        }

        for(int e = 0; e < nvir; ++e){
            for(int m = 0; m < nocc; ++m){
                for(int i = 0; i < nocc; ++i){
                    Fmi->add(m,i, (0.5)*t1->get(i,e)*fMat->get(m,e+nocc));
                }
            }
        }

        for(int e = 0; e < nvir; ++e){
            for(int n = 0; n < nocc; ++n){
                for(int m = 0; m < nocc; ++m){
                    for(int i = 0; i < nocc; ++i){
                        Fmi->add(m,i, spo_TEI[m][n][i][e+nocc]*t1->get(n,e));
                    }
                }
            }
        }

        for(int n = 0; n < nocc; ++n){
            for(int e = 0; e < nvir; ++e){
                for(int f = 0; f < nvir; ++f){
                    for(int m = 0; m < nocc; ++m){
                        for(int i = 0; i < nocc; ++i){
                            Fmi->add(m,i,(0.5)*taut[i][n][e][f]*spo_TEI[m][n][e+nocc][f+nocc]);
                        }
                    }
                }
            }
        }

        //Build Fme (nocc X nvir)

        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                Fme->set(m,e, fMat->get(m,e+nocc));
            }
        }
        for(int n = 0; n < nocc; ++n){
            for(int f = 0; f < nvir; ++f){
                for(int m = 0; m < nocc; ++m){
                    for(int e = 0; e < nvir; ++e){
                        Fme->add(m,e, t1->get(n,f)*spo_TEI[m][n][e+nocc][f+nocc]);
                    }
                }
            }
        }


        /*Now for the two-particle intermediates*/


        //Build Wmnij
        quad Wmnij(boost::extents[nocc][nocc][nocc][nocc]);

        for(int m = 0; m < nocc; ++m){
            for(int n = 0; n < nocc; ++n){
                for(int i = 0; i < nocc; ++i){
                    for(int j = 0; j < nocc; ++j){
                        Wmnij[m][n][i][j] = spo_TEI[m][n][i][j];
                    }
                }
            }
        }

        for(int e = 0; e < nvir; ++e){
            for(int m = 0; m < nocc; ++m){
                for(int n = 0; n < nocc; ++n){
                    for(int i = 0; i < nocc; ++i){
                        for(int j = 0; j < nocc; ++j){
                            Wmnij[m][n][i][j] += t1->get(j,e)*spo_TEI[m][n][i][e+nocc];
                        }
                    }
                }
            }
        }

        for(int e = 0; e < nvir; ++e){
            for(int m = 0; m < nocc; ++m){
                for(int n = 0; n < nocc; ++n){
                    for(int i = 0; i < nocc; ++i){
                        for(int j = 0; j < nocc; ++j){
                            Wmnij[m][n][i][j] += (-1.0)*t1->get(i,e)*spo_TEI[m][n][j][e+nocc];
                        }
                    }
                }
            }
        }

        for(int e = 0; e < nvir; ++e){
            for(int f = 0; f < nvir; ++f){
                for(int m = 0; m < nocc; ++m){
                    for(int n = 0; n < nocc; ++n){
                        for(int i = 0; i < nocc; ++i){
                            for(int j = 0; j < nocc; ++j){
                                Wmnij[m][n][i][j] += (0.25)*tau[i][j][e][f]*spo_TEI[m][n][e+nocc][f+nocc];
                            }
                        }
                    }
                }
            }
        }

        //Build Wabef (nvir X nvir X nvir X nvir )

        quad Wabef(boost::extents[nvir][nvir][nvir][nvir]);

        for(int a = 0; a < nvir; ++a){
            for(int b = 0; b < nvir; ++b){
                for(int e = 0; e < nvir; ++e){
                    for(int f = 0; f < nvir; ++f){
                        Wabef[a][b][e][f] = spo_TEI[a+nocc][b+nocc][e+nocc][f+nocc];
                    }
                }
            }
        }

        for(int m = 0; m < nocc; ++m){
            for(int a = 0; a < nvir; ++a){
                for(int b = 0; b < nvir; ++b){
                    for(int e = 0; e < nvir; ++e){
                        for(int f = 0; f < nvir; ++f){
                            Wabef[a][b][e][f] += (-1.0)*t1->get(m,b)*spo_TEI[a+nocc][m][e+nocc][f+nocc];
                        }
                    }
                }
            }
        }

        for(int m = 0; m < nocc; ++m){
            for(int a = 0; a < nvir; ++a){
                for(int b = 0; b < nvir; ++b){
                    for(int e = 0; e < nvir; ++e){
                        for(int f = 0; f < nvir; ++f){
                            Wabef[a][b][e][f] += t1->get(m,a)*spo_TEI[b+nocc][m][e+nocc][f+nocc];
                        }
                    }
                }
            }
        }

        for(int m = 0; m < nocc; ++m){
            for(int n = 0; n < nocc; ++n){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        for(int e = 0; e < nvir; ++e){
                            for(int f = 0; f < nvir; ++f){
                                Wabef[a][b][e][f] += (0.25)*tau[m][n][a][b]*spo_TEI[m][n][e+nocc][f+nocc];
                            }
                        }
                    }
                }
            }
        }

        //Build Wmbej (nocc x nvir x nvir x nocc)

        quad Wmbej(boost::extents[nocc][nvir][nvir][nocc]);

        for(int m = 0; m < nocc; ++m){
            for(int b = 0; b < nvir; ++b){
                for(int e = 0; e < nvir; ++e){
                    for(int j = 0; j < nocc; ++j){
                        Wmbej[m][b][e][j] = spo_TEI[m][b+nocc][e+nocc][j];
                    }
                }
            }
        }

        for(int f = 0; f < nvir; ++f){
            for(int m = 0; m < nocc; ++m){
                for(int b = 0; b < nvir; ++b){
                    for(int e = 0; e < nvir; ++e){
                        for(int j = 0; j < nocc; ++j){
                            Wmbej[m][b][e][j] += t1->get(j,f)*spo_TEI[m][b+nocc][e+nocc][f+nocc];
                        }
                    }
                }
            }
        }

        for(int n = 0; n < nocc; ++n){
            for(int b = 0; b < nvir; ++b){
                for(int m = 0; m < nocc; ++m){
                    for(int e = 0; e < nvir; ++e){
                        for(int j = 0; j < nocc; ++j){
                            Wmbej[m][b][e][j] += (-1.0)*t1->get(n,b)*spo_TEI[m][n][e+nocc][j];
                        }
                    }
                }
            }
        }

        for(int n = 0; n < nocc; ++n){
            for(int f = 0; f < nvir; ++f){
                for(int j = 0; j < nocc; ++j){
                    for(int b = 0; b < nvir; ++b){
                        for(int m = 0; m < nocc; ++m){
                            for(int e = 0; e < nvir; ++e){
                                Wmbej[m][b][e][j] += (-0.5)*t2[j][n][f][b]*spo_TEI[m][n][e+nocc][f+nocc] - t1->get(j,f)*t1->get(n,b)*spo_TEI[m][n][e+nocc][f+nocc];
                            }
                        }
                    }
                }
            }
        }

        /*Now that the F and W intermediates are defined, we need the updated t1 and t2 amplitudes*/

        //t1 amplitudes

        boost::shared_ptr<Matrix> t1D(new Matrix("t1D", nocc, nvir));
        t1D->zero();

        //first term
        for(int i = 0; i < nocc; ++i){
            for(int a = 0; a < nvir; ++a){
                t1D->set(i,a,fMat->get(i,a+nocc));
            }
        }

        //second term
        for(int e = 0; e < nvir; ++e){
            for(int a = 0; a < nvir; ++a){
                for(int i = 0; i < nocc; ++i){
                    t1D->add(i,a,t1->get(i,e)*Fae->get(a,e));
                }
            }
        }


        //third term
        for(int m = 0; m < nocc; ++m){
            for(int i = 0; i < nocc; ++i){
                for(int a = 0; a < nvir; ++a){
                    t1D->add(i, a, (-1.0)*t1->get(m,a)*Fmi->get(m,i));
                }
            }
        }

        //fourth term
        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int i = 0; i < nocc; ++i){
                    for(int a = 0; a < nvir; ++a){
                        t1D->add(i,a,t2[i][m][a][e]*Fme->get(m,e));
                    }
                }
            }
        }

        //fifth term
        for(int n = 0; n < nocc; ++n){
            for(int f = 0; f < nvir; ++f){
                for(int i = 0; i < nocc; ++i){
                    for(int a = 0; a < nvir; ++a){
                        t1D->add(i,a,(-1.0)*t1->get(n,f)*spo_TEI[n][a+nocc][i][f+nocc]);
                    }
                }
            }
        }

        //sixth term
        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int f = 0; f < nvir; ++f){
                    for(int i = 0; i < nocc; ++i){
                        for(int a = 0; a < nvir; ++a){
                            t1D->add(i,a, (-0.5)*t2[i][m][e][f]*spo_TEI[m][a+nocc][e+nocc][f+nocc]);
                        }
                    }
                }
            }
        }

        //seventh term
        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int n = 0; n < nocc; ++n){
                    for(int i = 0; i < nocc; ++i){
                        for(int a = 0; a < nvir; ++a){
                            t1D->add(i,a, (-0.5)*t2[m][n][a][e]*spo_TEI[n][m][e+nocc][i]);
                        }
                    }
                }
            }
        }


        //Now for the t2 equation

        //Define a 4D array (t2D) that will represent the t2 aplitudes multiplied by the two-particle
        //denominators.

        quad t2D(boost::extents[nocc][nocc][nvir][nvir]);

        //All matrices of the form t2D# are used as intermediate terms in the t2 equation. Hopefully
        //I'll reduce the number of these

        boost::shared_ptr<Matrix> t2D1(new Matrix("t2D1", nvir, nvir));

        t2D.empty();
        t2D1->zero();

        //first term
        for(int i = 0; i < nocc; ++i){
            for(int j = 0; j < nocc; ++j){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        t2D[i][j][a][b] = spo_TEI[i][j][a+nocc][b+nocc];
                    }
                }
            }
        }

        //second term
        boost::shared_ptr<Matrix> t2D2(new Matrix("t2D2", nvir, nvir));
        t2D2->zero();

        for(int m = 0; m < nocc; ++m){
            for(int b = 0; b < nvir; ++b){
                for(int e = 0; e < nvir; ++e){
                    t2D1->add(b,e, (0.5)*t1->get(m,b)*Fme->get(m,e));
                }
            }
        }

        for(int b = 0; b < nvir; ++b){
            for(int e = 0; e < nvir; ++e){
                t2D2->set(b,e, Fae->get(b,e) - t2D1->get(b,e));
            }
        }

        for(int e = 0; e < nvir; ++e){
            for(int i = 0; i < nocc; ++i){
                for(int j = 0; j < nocc; ++j){
                    for(int a = 0; a < nvir; ++a){
                        for(int b = 0; b < nvir; ++b){
                            t2D[i][j][a][b] += t2[i][j][a][e]*t2D2->get(b,e);
                        }
                    }
                }
            }
        }

        //third term
        boost::shared_ptr<Matrix> t2D3(new Matrix("t2D3", nvir, nvir));
        t2D1->zero();
        t2D3->zero();

        for(int m = 0; m < nocc; ++m){
            for(int b = 0; b < nvir; ++b){
                for(int e = 0; e < nvir; ++e){
                    t2D1->add(b,e, (0.5)*t1->get(m,b)*Fme->get(m,e));
                }
            }
        }

        for(int b = 0; b < nvir; ++b){
            for(int e = 0; e < nvir; ++e){
                t2D2->set(b,e, Fae->get(b,e) - t2D1->get(b,e));
            }
        }

        for(int e = 0; e < nvir; ++e){
            for(int i = 0; i < nocc; ++i){
                for(int j = 0; j < nocc; ++j){
                    for(int a = 0; a < nvir; ++a){
                        for(int b = 0; b < nvir; ++b){
                            t2D[i][j][a][b] -= t2[i][j][b][e]*t2D2->get(a,e);
                        }
                    }
                }
            }
        }

        //fourth term
        boost::shared_ptr<Matrix> t2D4(new Matrix("t2D4", nocc, nocc));
        boost::shared_ptr<Matrix> t2D5(new Matrix("t2D5", nocc, nocc));
        t2D4->zero();
        t2D5->zero();

        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int j = 0; j < nocc; ++j){
                    t2D4->add(m,j, (0.5)*t1->get(j,e)*Fme->get(m,e));
                }
            }
        }

        for(int m = 0; m < nocc; ++m){
            for(int j = 0; j < nocc; ++j){
                t2D5->set(m,j, Fmi->get(m,j) + t2D4->get(m,j));
            }
        }

        for(int m = 0; m < nocc; ++m){
            for(int i = 0; i < nocc; ++i){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        for(int j = 0; j < nocc; ++j){
                            t2D[i][j][a][b] -= t2[i][m][a][b]*t2D5->get(m,j);
                        }
                    }
                }
            }
        }

        //fifth term
        t2D4->zero();
        t2D5->zero();

        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int j = 0; j < nocc; ++j){
                    t2D4->add(m,j, (0.5)*t1->get(j,e)*Fme->get(m,e));
                }
            }
        }

        for(int m = 0; m < nocc; ++m){
            for(int j = 0; j < nocc; ++j){
                t2D5->set(m,j, Fmi->get(m,j) + t2D4->get(m,j));
            }
        }

        for(int m = 0; m < nocc; ++m){
            for(int i = 0; i < nocc; ++i){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        for(int j = 0; j < nocc; ++j){
                            t2D[i][j][a][b] += t2[j][m][a][b]*t2D5->get(m,i);
                        }
                    }
                }
            }
        }

        //sixth term
        for(int m = 0; m < nocc; ++m){
            for(int n = 0; n < nocc; ++n){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        for(int i = 0; i < nocc; ++i){
                            for(int j = 0; j < nocc; ++j){
                                t2D[i][j][a][b] += (0.5)*tau[m][n][a][b]*Wmnij[m][n][i][j];
                            }
                        }
                    }
                }
            }
        }

        //seventh term
        for(int e = 0; e < nvir; ++e){
            for(int f = 0; f < nvir; ++f){
                for(int i = 0; i < nocc; ++i){
                    for(int j = 0; j < nocc; ++j){
                        for(int a = 0; a < nvir; ++a){
                            for(int b = 0; b < nvir; ++b){
                                t2D[i][j][a][b] += (0.5)*tau[i][j][e][f]*Wabef[a][b][e][f];
                            }
                        }
                    }
                }
            }
        }

        //eigth term
        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int i = 0; i < nocc; ++i){
                    for(int a = 0; a < nvir; ++a){
                        for(int b = 0; b < nvir; ++b){
                            for(int j = 0; j < nocc; ++j){
                                t2D[i][j][a][b] += t2[i][m][a][e]*Wmbej[m][b][e][j] - t1->get(i,e)*t1->get(m,a)*spo_TEI[m][b+nocc][e+nocc][j];
                            }
                        }
                    }
                }
            }
        }

        //ninth term
        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int i = 0; i < nocc; ++i){
                    for(int a = 0; a < nvir; ++a){
                        for(int b = 0; b < nvir; ++b){
                            for(int j = 0; j < nocc; ++j){
                                t2D[i][j][a][b] -= t2[j][m][a][e]*Wmbej[m][b][e][i] - t1->get(j,e)*t1->get(m,a)*spo_TEI[m][b+nocc][e+nocc][i];
                            }
                        }
                    }
                }
            }
        }

        //tenth term
        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int i = 0; i < nocc; ++i){
                    for(int a = 0; a < nvir; ++a){
                        for(int b = 0; b < nvir; ++b){
                            for(int j = 0; j < nocc; ++j){
                                t2D[i][j][a][b] -= t2[i][m][b][e]*Wmbej[m][a][e][j] - t1->get(i,e)*t1->get(m,b)*spo_TEI[m][a+nocc][e+nocc][j];
                            }
                        }
                    }
                }
            }
        }

        //eleventh term
        for(int m = 0; m < nocc; ++m){
            for(int e = 0; e < nvir; ++e){
                for(int i = 0; i < nocc; ++i){
                    for(int a = 0; a < nvir; ++a){
                        for(int b = 0; b < nvir; ++b){
                            for(int j = 0; j < nocc; ++j){
                                t2D[i][j][a][b] += t2[j][m][b][e]*Wmbej[m][a][e][i] - t1->get(j,e)*t1->get(m,b)*spo_TEI[m][a+nocc][e+nocc][i];
                            }
                        }
                    }
                }
            }
        }

        //twelfth term

        for(int e = 0; e < nvir; ++e){
            for(int i = 0; i < nocc; ++i){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        for(int j = 0; j < nocc; ++j){
                            t2D[i][j][a][b] += t1->get(i,e)*spo_TEI[a+nocc][b+nocc][e+nocc][j];
                        }
                    }
                }
            }
        }

        //thirteenth term

        for(int e = 0; e < nvir; ++e){
            for(int i = 0; i < nocc; ++i){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        for(int j = 0; j < nocc; ++j){
                            t2D[i][j][a][b] -= t1->get(j,e)*spo_TEI[a+nocc][b+nocc][e+nocc][i];
                        }
                    }
                }
            }
        }

        //fourteenth term
        for(int m = 0; m < nocc; ++m){
            for(int a = 0; a < nvir; ++a){
                for(int b = 0; b < nvir; ++b){
                    for(int i = 0; i < nocc; ++i){
                        for(int j = 0; j < nocc; ++j){
                            t2D[i][j][a][b] -= t1->get(m,a)*spo_TEI[m][b+nocc][i][j];
                        }
                    }
                }
            }
        }

        //fifteenth (last) term
        for(int m = 0; m < nocc; ++m){
            for(int a = 0; a < nvir; ++a){
                for(int b = 0; b < nvir; ++b){
                    for(int i = 0; i < nocc; ++i){
                        for(int j = 0; j < nocc; ++j){
                            t2D[i][j][a][b] += t1->get(m,b)*spo_TEI[m][a+nocc][i][j];
                        }
                    }
                }
            }
        }

        //Construct 2-particle denominator array
        quad D2(boost::extents[nocc][nocc][nvir][nvir]);
        D2.empty();

        for(int i = 0; i < nocc; ++i){
            for(int j = 0; j < nocc; ++j){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        D2[i][j][a][b] = fMat->get(i,i) + fMat->get(j,j) - fMat->get(a+nocc,a+nocc) - fMat->get(b+nocc,b+nocc);
                    }
                }
            }
        }

        //Now we need to take our terms t1D and t2D, and divide by the respective denominators

        //The updated t1 amps are:

        for(int i = 0; i < nocc; ++i){
            for(int a = 0; a < nvir; ++a){
                t1->set(i,a, t1D->get(i,a) / (fMat->get(i,i) - fMat->get(a+nocc, a+nocc)) );
            }
        }

        //The updated t2 amps are:

        for(int i = 0; i < nocc; ++i){
            for(int j = 0; j < nocc; ++j){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        t2prev[i][j][a][b] = t2[i][j][a][b];
                        t2[i][j][a][b] = t2D[i][j][a][b]/D2[i][j][a][b];
                    }
                }
            }
        }

        //Transform amps using SRG equations
        if(options.get_bool("SRG") == true){
            for(int i = 0; i < nocc; ++i){
                for(int j = 0; j < nocc; ++j){
                    for(int a = 0; a < nvir; ++a){
                        for(int b = 0; b < nvir; ++b){
                            t2[i][j][a][b] = t2[i][j][a][b]*( 1.0 - exp(-1.0*D2[i][j][a][b]*D2[i][j][a][b]*options.get_double("srg_s")) );
                        }
                    }
                }
            }

            for(int i = 0; i < nocc; ++i){
                for(int a = 0; a < nvir; ++a){
                    t1->set(i,a, t1->get(i,a)*( 1.0 - exp(-1.0*options.get_double("srg_s")*(fMat->get(i,i) - fMat->get(a+nocc, a+nocc))*(fMat->get(i,i) - fMat->get(a+nocc, a+nocc)))));
                }
            }
        }

        //Compute the energy
        Eprev = Ecc;
        Ecc = 0.0;

        for(int i = 0; i < nocc; ++i){
            for(int a = 0; a < nvir; ++a){
                Ecc += fMat->get(i,a+nocc) * t1->get(i,a+nocc);
                for(int j = 0; j < nocc; ++j){
                    for(int b = 0; b < nvir; ++b){
                        Ecc += (0.25)*spo_TEI[i][j][a+nocc][b+nocc]*t2[i][j][a][b] + (0.5)*spo_TEI[i][j][a+nocc][b+nocc]*t1->get(i,a)*t1->get(j,b);
                    }
                }
            }
        }

        //Count the iterations
        count++;

        //Compute the energy difference wrt the previous iteration
        dE = fabs(Ecc - Eprev);

        //calculate RMS in t2 amps
        double RMS = 0.0;
        for(int i = 0; i < nocc; ++i){
            for(int j = 0; j < nocc; ++j){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        RMS += (t2[i][j][a][b] - t2prev[i][j][a][b])*(t2[i][j][a][b] - t2prev[i][j][a][b]);
                    }
                }
            }
        }
        RMS = sqrt(RMS);

        //Print the results for each iteration
        outfile->Printf("\n   [%3d]      %20.12f   %20.12f   %14.6f", count, Ecc, dE, RMS);

        //Convergence criteria
        if(dE < options.get_double("e_conv") and RMS < options.get_double("t2_conv") ){
            converged = true;
            outfile->Printf("\n\nCCSD energy has converged in %d cycles.\n", count);
        }

    }//end while loop

    double Et = 0.0;
    //Calclulate perturbative triples correction (CCSD(T))
    if(options.get_bool("P_TRIPLES") == true){
    outfile->Printf("\nComputing triples (T) correction\n");
    double Eint;

    for(int i = 0; i < nocc; ++i){
        for(int j = 0; j < nocc; ++j){
            for(int k = 0; k < nocc; ++k){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        for(int c = 0; c < nvir; ++c){
                            Eint = 0.0;
                            for(int e = 0; e < nvir; ++e){
                                Eint +=    (t2[j][k][a][e]*spo_TEI[e+nocc][i][b+nocc][c+nocc] - t2[j][k][b][e]*spo_TEI[e+nocc][i][a+nocc][c+nocc] - t2[j][k][c][e]*spo_TEI[e+nocc][i][b+nocc][a+nocc]
                                          - t2[i][k][a][e]*spo_TEI[e+nocc][j][b+nocc][c+nocc] + t2[i][k][b][e]*spo_TEI[e+nocc][j][a+nocc][c+nocc] + t2[i][k][c][e]*spo_TEI[e+nocc][j][b+nocc][a+nocc]
                                          - t2[j][i][a][e]*spo_TEI[e+nocc][k][b+nocc][c+nocc] + t2[j][i][b][e]*spo_TEI[e+nocc][k][a+nocc][c+nocc] + t2[j][i][c][e]*spo_TEI[e+nocc][k][b+nocc][a+nocc]);
                            }
                            for(int m = 0; m < nocc; ++m){
                                Eint += ((-1.0)*t2[i][m][b][c]*spo_TEI[m][a+nocc][j][k] + t2[i][m][a][c]*spo_TEI[m][b+nocc][j][k] + t2[i][m][b][a]*spo_TEI[m][c+nocc][j][k]
                                              + t2[j][m][b][c]*spo_TEI[m][a+nocc][i][k] - t2[j][m][a][c]*spo_TEI[m][b+nocc][i][k] - t2[j][m][b][a]*spo_TEI[m][c+nocc][i][k]
                                               + t2[k][m][b][c]*spo_TEI[m][a+nocc][j][i] - t2[k][m][a][c]*spo_TEI[m][b+nocc][j][i] - t2[k][m][b][a]*spo_TEI[m][c+nocc][j][i]);
                            }
                            Et +=  ((Eint*Eint) + ((t1->get(i,a)*spo_TEI[j][k][b+nocc][c+nocc] - t1->get(i,b)*spo_TEI[j][k][a+nocc][c+nocc] - t1->get(i,c)*spo_TEI[j][k][b+nocc][a+nocc]
                                   -t1->get(j,a)*spo_TEI[i][k][b+nocc][c+nocc] + t1->get(j,b)*spo_TEI[i][k][a+nocc][c+nocc] + t1->get(j,c)*spo_TEI[i][k][b+nocc][a+nocc]
                                   -t1->get(k,a)*spo_TEI[j][i][b+nocc][c+nocc] + t1->get(k,b)*spo_TEI[j][i][a+nocc][c+nocc] + t1->get(k,c)*spo_TEI[j][i][b+nocc][a+nocc] )*(Eint)))
                                    / (fMat->get(i,i) + fMat->get(j,j) + fMat->get(k,k) - fMat->get(a+nocc,a+nocc) - fMat->get(b+nocc,b+nocc) - fMat->get(c+nocc,c+nocc));
                        }
                    }
                }
            }
        }
    }

    Et /= 36.0;
    }//end if

    //Print a summary
    outfile->Printf("\nSCF energy                                  %20.15f a.u.", Eref);
    if(options.get_bool("SRG") == false){
        outfile->Printf("\nMP2 energy correction                       %20.15f a.u.", Emp2);
        outfile->Printf("\nMP2 total energy                            %20.15f a.u. \n\n", Eref+Emp2);
        outfile->Printf("CCSD correlation energy                     %20.15f a.u.", Ecc);
        outfile->Printf("\nCCSD total energy                           %20.15f a.u.", Ecc + Eref );
        if(options.get_bool("P_TRIPLES") == true){
            outfile->Printf("\nPerturbative Triples (T) correction         %20.15f a.u.", Et);
            outfile->Printf("\nCCSD(T) total energy                        %20.15f a.u.", Et + Ecc + Eref);
        }
    }

    if(options.get_bool("SRG") == true){
        outfile->Printf("\nMP2 energy correction                       %20.15f a.u.", Emp2);
        outfile->Printf("\nMP2 total energy                            %20.15f a.u. \n\n", Eref+Emp2);
        outfile->Printf("SRG-CCSD correlation energy                     %20.15f a.u.", Ecc);
        outfile->Printf("\nSRG-CCSD total energy                           %20.15f a.u.", Ecc + Eref );
        if(options.get_bool("P_TRIPLES") == true){
            outfile->Printf("\nSRG Perturbative Triples (T) correction         %20.15f a.u.", Et);
            outfile->Printf("\nSRG-CCSD(T) total energy                        %20.15f a.u.", Et + Ecc + Eref);
        }
    }



    outfile->Flush();

    return Success;
}

}} // End Namespaces
