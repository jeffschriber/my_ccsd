/*
 * @BEGIN LICENSE
 *
 * new_ccsd by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */
#include <math.h>
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"


// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

namespace psi{ namespace new_ccsd{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "NEW_CCSD" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);

        options.add_double("CC_CONVERGENCE", 1e-9);
        options.add_double("T2_CONVERGENCE", 1e-5);

    }

    return true;
}

extern "C"
SharedWavefunction new_ccsd(SharedWavefunction ref_wfn, Options& options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");

    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Have the reference (SCF) wavefunction, ref_wfn
    if(!ref_wfn) throw PSIEXCEPTION("SCF has not been run yet!");

    // Quickly check that there are no open shell orbitals here...
    int nirrep  = ref_wfn->nirrep();

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example in the test suite.
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(ref_wfn, spaces, IntegralTransform::Restricted);
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

    // Let's store the antisymmetrized two-electron integrals in a std::vector for optimal speed
    int nmo = ref_wfn->nmo(); 
    int nmo2 = nmo * nmo;
    int nmo3 = nmo2 * nmo;
    int nmo4 = nmo3 * nmo;
    std::vector<double> spat_tei(nmo4); 


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
//                psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f, "
//                                 "symmetries = (%1d %1d | %1d %1d), "
//                                 "relative indices = (%2d %2d | %2d %2d)\n",
//                                 p, q, r, s, K.matrix[h][pq][rs],
//                                 psym, qsym, rsym, ssym,
//                                 prel, qrel, rrel, srel);

                // Now fill the tei vector 
                spat_tei[nmo3*p + nmo2*q + nmo*r + s] = K.matrix[h][pq][rs];
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);


    // Define number of spin orbitals
    int nso = nmo*2;
    int nso2 = nso*nso;
    int nso3 = nso2*nso;
    int nso4 = nso3*nso;

    // Another vector for the teis in the spin orbital basis
    std::vector<double> tei(nso4, 0.0);

    // Lets define a quick function to easily index the integrals
    auto four_idx = [&](int p, int q, int r, int s) -> int 
    {
        return (p*nso3 + q*nso2 + r*nso + s); 
    };

    // Now lets build the antisymmetrized teis
    // We'll alternate between alpha and beta
    for( int p = 0; p < nso; ++p){
        for( int q = 0; q < nso; ++q){
            for( int r = 0; r < nso; ++r){
                for( int s = 0; s < nso; ++s){
                    int pp = p / 2; 
                    int qq = q / 2; 
                    int rr = r / 2; 
                    int ss = s / 2; 
    
                    double v1 = spat_tei[nmo3*pp + nmo2*rr + nmo*qq + ss] * (p%2==r%2) * (q%2==s%2); 
                    double v2 = spat_tei[nmo3*pp + nmo2*ss + nmo*qq + rr] * (p%2==s%2) * (q%2==r%2);
                    tei[four_idx(p,q,r,s)] = v1 - v2;
                }
            }
        }
    }

    // Now get the AO fock Matrix and transform it
    SharedMatrix Fao = ref_wfn->Fa();
    SharedMatrix Ca = ref_wfn->Ca();
    
    //This will be the Fock matrix in the MO basis
    SharedMatrix Fmo(new Matrix(nmo, nmo));
    Fmo->transform( Ca, Fao, Ca);

    //This will be the Fock matrix in the SO basis
    std::vector<double> F(nso);
    for( int i = 0; i < nmo; ++i ){
        F[i*2] = Fmo->get(i,i);
        F[i*2 + 1] = Fmo->get(i,i);
    }

    // Get dimensions of occ and vir
    int nocc = 2 * ref_wfn->nalpha();
    int nvir = nso - nocc;


    // Build initial t1 and t2 amplitudes
    double t1[nocc][nvir];
    double t2[nocc][nocc][nvir][nvir];


    for( int i = 0; i < nocc; ++i ){
        for( int a = 0; a < nvir; ++a ){
            t1[i][a] = 0.0;
        }
    }

    for( int i = 0; i < nocc; ++i){
        for( int j = 0; j < nocc; ++j){
            for( int a = nocc; a < nso; ++a){
                for( int b = nocc; b < nso; ++b){
                    t2[i][j][a-nocc][b-nocc] = tei[four_idx(i,j,a,b)] / (F[i] + F[j] - F[a] - F[b]);
                }
            }
        }
    }


    // Compute and print mp2 energy
    double Emp2 = 0.0;

    for( int i = 0; i < nocc; ++i){
        for( int j = 0; j < nocc; ++j){
            for( int a = nocc; a < nso; ++a){
                for( int b = nocc; b < nso; ++b){
                    Emp2 += tei[four_idx(i,j,a,b)] * t2[i][j][a-nocc][b-nocc]; 
                }
            }
        }
    }
    Emp2 *= 0.25;

    // Begin CCSD procedure
    bool converged = false;
    double old_energy = 0.0;
    int count =  0;
    double Ecc = 0.0;

    outfile->Printf("\n\n  -------------");
    outfile->Printf(  "\n   CCSD Energy");
    outfile->Printf(  "\n     By Jeff");
    outfile->Printf("\n  -------------\n\n");
    

    outfile->Printf("  -----------------------------------------------------------------------------");
    outfile->Printf("\n      Iter                E(CCSD)                Delta(E)            RMS(T2)   \n");
    outfile->Printf("  -----------------------------------------------------------------------------");
    
    while ( !converged ){
        
        // The two-particle excitation operators
        double taut[nocc][nocc][nvir][nvir];
        double tau[nocc][nocc][nvir][nvir];
        for( int i = 0; i < nocc; ++i){
            for( int j = 0; j < nocc; ++j){
                for( int a = nocc; a < nso; ++a){
                    for( int b = nocc; b < nso; ++b){
                        taut[i][j][a-nocc][b-nocc] = t2[i][j][a-nocc][b-nocc] 
                                                     + 0.5*(t1[i][a-nocc]*t1[j][b-nocc] 
                                                     - t1[i][b-nocc]*t1[j][a-nocc]);
                        tau[i][j][a-nocc][b-nocc] = t2[i][j][a-nocc][b-nocc] 
                                                    + t1[i][a-nocc]*t1[j][b-nocc] 
                                                    - t1[i][b-nocc]*t1[j][a-nocc];
                    }
                }
            }
        }

        // The vir-vir F intermediate
        SharedMatrix Fvv(new Matrix(nvir,nvir));
double sumsq = 0.0;
        double value;
        for( int a = nocc; a < nso; ++a){
            for( int e = nocc; e < nso; ++e){
                value = 0.0;
                for( int m = 0; m < nocc; ++m ){
                    for( int f = nocc; f < nso; ++f ){
                        value += t1[m][f-nocc]*tei[four_idx(m,a,f,e)];
                    }
                }
                
                for( int m = 0; m < nocc; ++m ){
                    for( int n = 0; n < nocc; ++n ){
                        for( int f = nocc; f < nso; ++f ){
                            value -= 0.5 * taut[m][n][a-nocc][f-nocc]*tei[four_idx(m,n,e,f)]; 
                        }
                    }
                }
                Fvv->set(a-nocc,e-nocc,value);
                sumsq += value * value;
            }        
        }        
outfile->Printf("\n F_ae = %1.12f", sumsq);
        // The occ-occ F intermediate
sumsq=0.0;
        SharedMatrix Foo(new Matrix(nocc,nocc));
        for(int m = 0; m < nocc; ++m){
            for( int i = 0; i < nocc; ++i){
                value = 0.0;
                for( int e = nocc; e < nso; ++e ){
                    for( int n = 0; n < nocc; ++n ){
                        value += t1[n][e-nocc]*tei[four_idx(m,n,i,e)];
                    }
                }
                for( int e = nocc; e < nso; ++e ){
                    for( int f = nocc; f < nso; ++f ){
                        for( int n = 0; n < nocc; ++n ){
                            value += 0.5 * taut[i][n][e-nocc][f-nocc] * tei[four_idx(m,n,e,f)];
                        }
                    }
                }
                Foo->set(m,i,value); 
                sumsq += value * value;
            }
        } 
outfile->Printf("\n F_mi = %1.12f", sumsq);

        // The occ-vir F intermediate
sumsq=0.0;
        SharedMatrix Fov(new Matrix(nocc, nvir));
        for( int m = 0; m < nocc; ++m ){
            for( int e = nocc; e < nso; ++e){
                value = 0.0;
                for( int n = 0; n < nocc; ++n ){
                    for( int f = nocc; f < nso; ++f ){
                        value += t1[n][f-nocc]*tei[four_idx(m,n,e,f)];    
                    }
                }
                Fov->set(m,e-nocc, value);
                sumsq += value * value;
            }
        }
outfile->Printf("\n F_me = %1.12f", sumsq);

        // The occ-occ-occ-occ W intermediate
    double sumsq1 = 0.0;
        double Woooo[nocc][nocc][nocc][nocc];
        
        for( int m = 0; m < nocc; ++m ){
            for( int n = 0; n < nocc; ++n ){
                for( int i = 0; i < nocc; ++i ){
                    for( int j = 0; j < nocc; ++j ){
                        value = tei[four_idx(m,n,i,j)]; 
                        
                        for( int e = nocc; e < nso; ++e ){
                            value += t1[j][e-nocc] * tei[four_idx(m,n,i,e)]
                                   - t1[i][e-nocc] * tei[four_idx(m,n,j,e)];
                        }

                        for( int e = nocc; e < nso; ++e ){
                            for( int f = nocc; f < nso; ++f ){
                                value += 0.25 * tau[i][j][e-nocc][f-nocc] * tei[four_idx(m,n,e,f)]; 
                            }   
                        }
                        Woooo[m][n][i][j] = value;
                        sumsq1 += (value*value);
                    }    
                }    
            }    
        }    
outfile->Printf("\n  Wmnij: %1.9f", sumsq1);
        // The vir-vir-vir-vir W intermediate
    double sumsq2 = 0.0;
        double Wvvvv[nvir][nvir][nvir][nvir];
        for(int a = nocc; a < nso; ++a ){
            for(int b = nocc; b < nso; ++b ){
                for(int e = nocc; e < nso; ++e ){
                    for(int f = nocc; f < nso; ++f ){
                        value = tei[four_idx(a,b,e,f)];

                        for( int m = 0; m < nocc; ++m ){
                            value -= (t1[m][b-nocc] * tei[four_idx(a,m,e,f)]
                                    - t1[m][a-nocc] * tei[four_idx(b,m,e,f)]);
                        }
                        for( int m = 0; m < nocc; ++m ){
                            for( int n = 0; n < nocc; ++n ){
                                value += 0.25 * tau[m][n][a-nocc][b-nocc] * tei[four_idx(m,n,e,f)];
                            }
                        }
                        Wvvvv[a-nocc][b-nocc][e-nocc][f-nocc] = value;
                        sumsq2 += (value*value);
                    }
                }
            }
        }
outfile->Printf("\n  Wabef: %1.9f", sumsq2);

        // The occ-vir-vir-occ W intermediate
    double sumsq3 = 0.0;
        double Wovvo[nocc][nvir][nvir][nocc];
        for( int m = 0; m < nocc; ++m ){
            for( int b = nocc; b < nso; ++b ){
                for( int e = nocc; e < nso; ++e ){
                    for( int j = 0; j < nocc; ++j ){
                        value = tei[four_idx(m,b,e,j)];

                        for( int f = nocc; f < nso; ++f ){
                            value += t1[j][f-nocc] * tei[four_idx(m,b,e,f)];
                        }
                        for( int n = 0; n < nocc; ++n ){
                            value -= t1[n][b-nocc] * tei[four_idx(m,n,e,j)];
                        }
                        for( int f = nocc; f < nso; ++f ){
                            for( int n = 0; n < nocc; ++n ){
                                value -= (0.5*t2[j][n][f-nocc][b-nocc] + t1[j][f-nocc]*t1[n][b-nocc]) *
                                          tei[four_idx(m,n,e,f)];
                            }
                        }
                        Wovvo[m][b-nocc][e-nocc][j] = value;
                        sumsq3 += (value*value);
                    }
                }
            }
        }
outfile->Printf("\n  Wabej: %1.9f", sumsq3);

        // Compute updated T1 amplitudes
        double t1new[nocc][nvir];
double t1sq = 0.0;
        for( int i = 0; i < nocc; ++i ){
            for( int a = nocc; a < nso; ++a ){
        value = 0.0;
                
                for( int e = nocc; e < nso; ++e ){
                    value += t1[i][e-nocc] * Fvv->get(a-nocc,e-nocc); 
                }
                for( int m = 0; m < nocc; ++m ){
                    value -= t1[m][a-nocc] * Foo->get(m,i);
                }
                for( int m = 0; m < nocc; ++m ){
                    for( int e = nocc; e < nso; ++e ){
                        value += t2[i][m][a-nocc][e-nocc] * Fov->get(m,e-nocc);
                    }
                }
                for( int n = 0; n < nocc; ++n ){
                    for( int f = nocc; f < nso; ++f ){    
                        value -= t1[n][f-nocc] * tei[four_idx(n,a,i,f)];
                    }
                }
                for( int m = 0; m < nocc; ++m ){
                    for( int e = nocc; e < nso; ++e ){
                        for( int f = nocc; f < nso; ++f ){
                            value -= 0.5 * t2[i][m][e-nocc][f-nocc] * tei[four_idx(m,a,e,f)];   
                        }
                    }
                }
                for( int m = 0; m < nocc; ++m ){
                    for( int e = nocc; e < nso; ++e ){
                        for( int n = 0; n < nocc; ++n ){
                            value -= 0.5 * t2[m][n][a-nocc][e-nocc] * tei[four_idx(n,m,e,i)];   
                        }
                    }
                }
t1sq += (value / (F[i] - F[a])) * (value / (F[i] - F[a]));
                t1new[i][a-nocc] = value / (F[i] - F[a]);
//outfile->Printf("\n t1[%d][%d] = %1.8f", i, a-nocc, value / (F[i] - F[a]));
            }
        }
outfile->Printf("\n  T1: %1.12f", t1sq);

        // Compute updated T2 amplitudes
        double t2new[nocc][nocc][nvir][nvir];
double t2sq = 0.0;
        for( int i = 0; i < nocc; ++i ){
            for( int j = 0; j < nocc; ++j ){
                for( int a = nocc; a < nso; ++a ){
                    for( int b = nocc; b < nso; ++b ){
                        value = tei[four_idx(i,j,a,b)];   
                        for( int e = nocc; e < nso; ++e ){
                            double tmp = 0.0;
                            for( int m = 0; m < nocc; ++m ){
                                tmp += t1[m][b-nocc] * Fov->get(m,e-nocc);
                            }
                            tmp *= 0.5;
                            value += t2[i][j][a-nocc][e-nocc] * (Fvv->get(b-nocc,e-nocc) - tmp);
                        }
                        for( int e = nocc; e < nso; ++e ){
                            double tmp = 0.0;
                            for( int m = 0; m < nocc; ++m ){
                                tmp += t1[m][a-nocc] * Fov->get(m,e-nocc);
                            }
                            tmp *= 0.5;
                            value -= t2[i][j][b-nocc][e-nocc] * (Fvv->get(a-nocc,e-nocc) - tmp);
                        }
                        
                        for( int m = 0; m < nocc; ++m ){
                            double tmp = 0.0;
                            for( int e = nocc; e < nso; ++e ){
                                tmp += t1[j][e-nocc] * Fov->get(m,e-nocc);
                            }
                            tmp *= 0.5;
                            value -= t2[i][m][a-nocc][b-nocc] * (Foo->get(m,j) + tmp);
                        }

                        for( int m = 0; m < nocc; ++m ){
                            double tmp = 0.0;
                            for( int e = nocc; e < nso; ++e ){
                                tmp += t1[i][e-nocc] * Fov->get(m,e-nocc);
                            }
                            tmp *= 0.5;
                            value += t2[j][m][a-nocc][b-nocc] * (Foo->get(m,i) + tmp);
                        }
                        
                        for( int m = 0; m < nocc; ++m ){
                            for( int n = 0; n < nocc; ++n ){
                                value += 0.5 * tau[m][n][a-nocc][b-nocc] * Woooo[m][n][i][j];
                            }
                        } 
                        for( int e = 0; e < nvir; ++e ){
                            for( int f = 0; f < nvir; ++f ){
                                value += 0.5 * tau[i][j][e][f] * Wvvvv[a-nocc][b-nocc][e][f];
                            }
                        } 

                        for( int m = 0; m < nocc; ++m ){
                            for( int e = nocc; e < nso; ++e ){
                                value += t2[i][m][a-nocc][e-nocc] * Wovvo[m][b-nocc][e-nocc][j] 
                                        - t1[i][e-nocc]*t1[m][a-nocc]*tei[four_idx(m,b,e,j)]; 
                                value -= t2[j][m][a-nocc][e-nocc] * Wovvo[m][b-nocc][e-nocc][i] 
                                        - t1[j][e-nocc]*t1[m][a-nocc]*tei[four_idx(m,b,e,i)]; 
                                value -= t2[i][m][b-nocc][e-nocc] * Wovvo[m][a-nocc][e-nocc][j] 
                                        - t1[i][e-nocc]*t1[m][b-nocc]*tei[four_idx(m,a,e,j)]; 
                                value += t2[j][m][b-nocc][e-nocc] * Wovvo[m][a-nocc][e-nocc][i] 
                                        - t1[j][e-nocc]*t1[m][b-nocc]*tei[four_idx(m,a,e,i)]; 
                            }
                        }

                        for( int e = nocc; e < nso; ++e ){
                            value += t1[i][e-nocc]*tei[four_idx(a,b,e,j)];
                            value -= t1[j][e-nocc]*tei[four_idx(a,b,e,i)];
                        }
                        for( int m = 0; m < nocc; ++m ){
                            value -= t1[m][a-nocc]*tei[four_idx(m,b,i,j)];
                            value += t1[m][b-nocc]*tei[four_idx(m,a,i,j)];
                        }

                        value /= ( F[i] + F[j] - F[a] - F[b] );

                        t2new[i][j][a-nocc][b-nocc] = value;
                        t2sq += value * value;
                    }
                }
            }
        }
outfile->Printf("\n  T2: %1.12f", t2sq);
        Ecc = 0.0;
        // Compute the energy with the new amplitudes
        for ( int i = 0; i < nocc; ++i ){
            for ( int j = 0; j < nocc; ++j ){
                for( int a = nocc; a < nso; ++a ){
                    for( int b = nocc; b < nso; ++b ){
                        double ijab = tei[four_idx(i,j,a,b)];
                        Ecc += 0.25 * ijab * t2new[i][j][a-nocc][b-nocc];
                        Ecc += 0.5 * ijab * t1new[i][a-nocc] * t1new[j][b-nocc];
                    }
                }
            }
        }
        
        double dE = std::abs( Ecc - old_energy );

        //calculate RMS in t2 amps
        double RMS = 0.0;
        for(int i = 0; i < nocc; ++i){
            for(int j = 0; j < nocc; ++j){
                for(int a = 0; a < nvir; ++a){
                    for(int b = 0; b < nvir; ++b){
                        RMS += (t2[i][j][a][b] - t2new[i][j][a][b])*(t2[i][j][a][b] - t2new[i][j][a][b]);
                    }
                }
            }
        }
        RMS = sqrt(RMS);

        //Print the results for each iteration
        outfile->Printf("\n     [%3d]      %20.12f   %20.12f   %14.6f", count, Ecc, dE, RMS);

        //Convergence criteria
        if(dE < options.get_double("CC_CONVERGENCE") and RMS < options.get_double("T2_CONVERGENCE") ){
            converged = true;
            outfile->Printf("\n\n  CCSD energy has converged in %d cycles.", count);
        }

        // Update the cluster amplitudes
        count++;
        old_energy = Ecc;
        std::copy(&t1new[0][0], &t1new[0][0] + nocc*nvir, &t1[0][0] );
        std::copy(&t2new[0][0][0][0], &t2new[0][0][0][0] + nocc*nocc*nvir*nvir, &t2[0][0][0][0] );
    } // End iterations

    // Print Summary
    double E_rhf = ref_wfn->reference_energy();

    outfile->Printf("\n\n  Energy Summary");
    outfile->Printf("\n  -----------------------------------------------------------");  
    outfile->Printf("\n     RHF energy (Eh):              %23.12f", E_rhf); 
    outfile->Printf("\n     MP2 correlation energy (Eh):  %23.12f", Emp2); 
    outfile->Printf("\n     MP2 total energy (Eh):        %23.12f", Emp2 + E_rhf); 
    outfile->Printf("\n     CCSD correlation energy (Eh): %23.12f", Ecc); 
    outfile->Printf("\n     CCSD total energy (Eh):       %23.12f", Ecc + E_rhf); 
    outfile->Printf("\n  -----------------------------------------------------------\n\n");  

    return ref_wfn;
}

}} // End Namespaces
