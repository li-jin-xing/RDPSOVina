/*
        PSOVina version 2.0

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. SIU <shirleysiu@umac.mo>


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

#include "rdpso_search.h"
#include "coords.h"
#include "quasi_newton.h"
#include "rdpso_mutate.h"
#include "mutate.h"

output_type pso_search::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_birds, double c1, double c2) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator, num_birds, c1, c2); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

// out is sorted
void pso_search::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_birds,double c1,double c2) const {
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);
	tmp.c.randomize(corner1, corner2, generator);  //first randomize
	fl best_e = max_fl;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	output_type tmp_rough = tmp;
	pso particle(num_birds, c1, c2,corner1,corner2,generator,tmp.c);
	double energy=0;
	int count=0, phrase=1,index[particle.number]={0};

	VINA_U_FOR(step, num_steps) {
		if(increment_me)
			++(*increment_me);
		output_type candidate = tmp_rough;
		output_type candidate_1 = tmp;

    if (phrase==1)
        rdpso_mutate_conf1(candidate, candidate_1, m, mutation_amplitude, generator, &particle, p, ig, g, hunt_cap, quasi_newton_par, step, num_steps, count);
    else
        rdpso_mutate_conf2(candidate, candidate_1, m, mutation_amplitude, generator, &particle, p, ig, g, hunt_cap, quasi_newton_par, step, num_steps, index);

    tmp_rough = candidate;

		if(candidate_1.e<0 && metropolis_accept(tmp.e, candidate_1.e, temperature, generator))
		{
  			tmp = candidate_1;
  			m.set(tmp.c); // FIXME? useless?
  			if(tmp.e < best_e || out.size() < num_saved_mins) {
  				quasi_newton_par(m, p, ig, tmp, g, authentic_v);
  				m.set(tmp.c); // FIXME? useless?
  				tmp.coords = m.get_heavy_atom_movable_coords();
  				add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
  				if(tmp.e < best_e)
  						best_e = tmp.e;
  			}
		}

		/***Criteria defined by RDPSOVina***/
		if(std::abs(pso::gbest_fit - energy) < 0.0001)
		{
  			count += 1;
  			if(phrase==1 && count>100)
  			{
				particle.sortParticle(index);
				phrase=2;
			}
			if(phrase==2 && count<3)
			{
				phrase=1;
				count=0;
			}
		}
    else
    {
  			energy = pso::gbest_fit;
  			count =0;
		}

	if(step>int(num_steps*0.9))
	{
		phrase=2;
	}
	}
	
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
