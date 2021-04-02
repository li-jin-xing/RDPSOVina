/*
       PSOVina version 2.0  

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

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

#include "rdpso_mutate.h"
#include <math.h>
#include <time.h>
#define PI 3.14159265

sz pso_count_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}

void rdpso_mutate_conf1(output_type& candidate, output_type& candidate_1, const model& m, fl amplitude, rng& generator, pso* particle, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par,int step, int num_steps, int count) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion

	output_type tmp_1 = candidate;
	output_type tmp_2 = candidate;

	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	srand((unsigned)time(NULL));
	double alpha_cm, beta_cm, isupdate;
	int y, mc_par=0;

	alpha_cm  = (1 - 0.3) * (1-step/(num_steps*0.9));
	if(alpha_cm < 0)
		alpha_cm=0;

	beta_cm  = 1.5 -  0.3 * (step/(num_steps*0.9));
	if(beta_cm < 1)
		beta_cm = 1;

	VINA_FOR_IN(i, candidate.c.ligands) {

		model tmp_m = m;
		const vec authentic_v(1000, 1000, 1000);

		mc_par = random_int(0, particle->number-1, generator);

		for (y=0;y<particle->number;y++)
		{
			candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
			candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
				candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);

			isupdate  = (double) rand() / (RAND_MAX + 1.0);
			if(isupdate<0.06)
			{
				quasi_newton_par(tmp_m, p, ig, candidate, g, v);
				particle->updateCurrentPosition(y,candidate.c.ligands[i].rigid.position);
				particle->updateCurrentOrientation(y,candidate.c.ligands[i].rigid.orientation);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					particle->updateCurrentTorsion(y, candidate.c.ligands[i].torsions[z],z);
			}
			else
			{
				if(particle->getPersonalBest(y)>100)
				{
					quasi_newton_par(tmp_m, p, ig, candidate, g, v);
					particle->updateCurrentPosition(y,candidate.c.ligands[i].rigid.position);
					particle->updateCurrentOrientation(y,candidate.c.ligands[i].rigid.orientation);
					for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
						particle->updateCurrentTorsion(y, candidate.c.ligands[i].torsions[z],z);
				}
				else
					candidate.e=tmp_m.eval_deriv(p, ig, v, candidate.c, g);
			}

			particle->updateCurrentBest(y,candidate.e);
			tmp_2 = candidate;


			if(y==mc_par && count<particle->number)
			{
				candidate.c.ligands[i].rigid.position =  particle->getPBestPosition(y);
				candidate.c.ligands[i].rigid.orientation = particle->getPBestOrientation(y);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate.c.ligands[i].torsions[z] = particle->getPBestTorsion(y,z);
				candidate.e = particle->getPersonalBest(y);
				markov_mutate_conf(candidate, tmp_m, 2, generator, p, ig, g, v, quasi_newton_par);
				if(candidate.e<tmp_2.e)
					tmp_2 = candidate;
				candidate_1=tmp_2;
			}


			if (tmp_2.e < particle->getPersonalBest(y))
			{
				particle->updatePersonalBest(y,tmp_2.e);
				particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
				particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
				for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
					particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
			}

			if(tmp_2.e < pso::gbest_fit)
			{
				particle->updateGlobalBest_1(tmp_2.e);
				pso::gbest_position = tmp_2.c.ligands[i].rigid.position;
				pso::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
				for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
					pso::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
			}
			if(count<50)
			{
				particle->computeRDPSOPositions(y,alpha_cm,beta_cm);
				fl gr = m.gyration_radius(i); 
				if(gr > epsilon_fl)
					particle->computeRDPSOOrientation(y,alpha_cm,beta_cm);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					particle->computeRDPSOTorsion(y,z,alpha_cm,beta_cm);
			}
			else
			{
				sz which1 = sz(which_int);
				VINA_CHECK(which1 < mutable_entities_num);
				
				if(which1 == 0){
					particle->computeRDPSOPositions(y,alpha_cm,beta_cm);
				}
				--which1;

				if(which1 == 0)
				{
					fl gr = m.gyration_radius(i); 
					if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
						particle->computeRDPSOOrientation(y,alpha_cm,beta_cm);
					}
				}
				--which1;

				if(which1 < candidate.c.ligands[i].torsions.size() && which1>=0)
				{
					particle->computeRDPSOTorsion(y,which1,alpha_cm,beta_cm);
				}
				
				which1 -= candidate.c.ligands[i].torsions.size();
			}
		}
	return;
	} 

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}

void rdpso_mutate_conf2(output_type& candidate, output_type& candidate_1, const model& m, fl amplitude, rng& generator, pso* particle, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par, int step, int num_steps, int* index) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion

	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	srand((unsigned)time(NULL));
	int y, par=random_int(0, 2, generator);
	VINA_FOR_IN(i, candidate.c.ligands) {

		model tmp_m = m;
		const vec authentic_v(1000, 1000, 1000);
		for (int x=0;x<3;x++)
		{
			y=index[x];
			candidate.c.ligands[i].rigid.position =  particle->getPBestPosition(y);
			candidate.c.ligands[i].rigid.orientation = particle->getPBestOrientation(y);
			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
				candidate.c.ligands[i].torsions[z] = particle->getPBestTorsion(y,z);
			markov_mutate_conf(candidate, tmp_m, 2, generator, p, ig, g, v, quasi_newton_par);
			if (candidate.e < particle->getPersonalBest(y))
			{
				particle->updatePersonalBest(y,candidate.e);
				particle->updateBestPosition(y,candidate.c.ligands[i].rigid.position);
				particle->updateBestOrientation(y,candidate.c.ligands[i].rigid.orientation);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					particle->updateBestTorsion(y, candidate.c.ligands[i].torsions[z],z);
			}
			if(candidate.e < pso::gbest_fit)
			{
				particle->updateGlobalBest_1(candidate.e);
				pso::gbest_position = candidate.c.ligands[i].rigid.position;
				pso::gbest_orientation = candidate.c.ligands[i].rigid.orientation;
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					pso::gbest_torsion[z] = candidate.c.ligands[i].torsions[z];
			}
			if(y==par && step>0.9*num_steps)
				candidate_1=candidate;
		}

	if(step>num_steps-10)
	{
		for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
			candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];
		candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;
		candidate_1.c.ligands[i].rigid.position = pso::gbest_position;
		candidate_1.e = pso::gbest_fit;
	}

	return; 
	} 

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}


void markov_mutate_conf(output_type& candidate, const model& m, fl amplitude, rng& generator, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par)
{ // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	output_type tmp_1 = candidate;
	output_type tmp_2 = candidate;
	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	int which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	VINA_FOR_IN(i, candidate.c.ligands)
	{
		model tmp_m = m;
		if(which == 0){
			candidate.c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator);
		}
		--which;

		if(which == 0)
		{
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
				vec rotation; 
				rotation = amplitude / gr * random_inside_sphere(generator); 
				quaternion_increment(candidate.c.ligands[i].rigid.orientation, rotation);
			}
		}
		--which;

		if(which < candidate.c.ligands[i].torsions.size() && which>=0) 
			candidate.c.ligands[i].torsions[which] = random_fl(-pi, pi, generator);
		
		which -= candidate.c.ligands[i].torsions.size();

		candidate.e=tmp_m.eval_deriv(p, ig, v, candidate.c, g);
		tmp_1 = candidate;
		quasi_newton_par(tmp_m, p, ig, tmp_1, g, v, 1);
		if(tmp_1.e < candidate.e)
		{
			candidate = tmp_1;
			tmp_2=candidate;
			quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
			if(tmp_2.e < candidate.e)
				candidate = tmp_2;
		}
	}

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}
