/*
        PSOVina version 1.0                     Date: 26/11/2014

        This file is revised from mutate.h in AutoDock Vina.

        Authors: Marcus C. K. Ng  <marcus.ckng@gmail.com>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_PSO_MUTATE_H
#define VINA_PSO_MUTATE_H

#include "rdpso.h"
#include "model.h"
#include "quasi_newton.h"

// does not set model
void rdpso_mutate_conf1(output_type& c, output_type& c_1, const model& m, fl amplitude, rng& generator, pso*, const precalculate&, const igrid&, change&, const vec&, quasi_newton&, int, int, int);
void rdpso_mutate_conf2(output_type& c, output_type& c_1, const model& m, fl amplitude, rng& generator, pso*, const precalculate&, const igrid&, change&, const vec&, quasi_newton&, int, int, int*);
void markov_mutate_conf(output_type& c, const model& m, fl amplitude, rng& generator, const precalculate&, const igrid&, change&, const vec&, quasi_newton&);

#endif
