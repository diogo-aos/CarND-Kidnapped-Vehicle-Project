/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define DEBUG 0
#define DEBUG_2 1


using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// values near the mean are the most likely
	// standard deviation affects the dispersion of generated values from the mean
	default_random_engine gen;
	std::normal_distribution<double> d_x(x, std[0]);
	std::normal_distribution<double> d_y(y, std[1]);
	std::normal_distribution<double> d_theta(theta, std[2]);

	num_particles = 100;
	for (int i=0; i<num_particles; i++){
		weights.push_back(1.0);
		Particle p;
		p.id = i;
		p.x = d_x(gen);
		p.y = d_y(gen);
		p.theta = d_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	std::normal_distribution<double> d_velocity(velocity, std_pos[0]);
	std::normal_distribution<double> d_yaw_rate(yaw_rate, std_pos[1]);

	std::normal_distribution<double> d_x(0, std_pos[0]);
	std::normal_distribution<double> d_y(0, std_pos[1]);
	std::normal_distribution<double> d_theta(0, std_pos[2]);

	for(int i=0; i<num_particles; i++){
		double x = particles[i].x,
				   y = particles[i].y,
					 theta = particles[i].theta;

		// double noisy_v = velocity + d_velocity(gen);
		// double noisy_yaw_rate = yaw_rate + d_yaw_rate(gen);

		double noisy_v = velocity;
		double noisy_yaw_rate = yaw_rate;

		if (fabs(yaw_rate) < 0.00001){  // yaw rate is 0
			x += noisy_v * delta_t * cos(theta);
			y += noisy_v * delta_t * sin(theta);
			// theta = theta;
		}
		else {  // yaw rate different than zero
			x += (noisy_v / noisy_yaw_rate) * ( sin(theta + noisy_yaw_rate * delta_t) - sin(theta) );
			y += (noisy_v / noisy_yaw_rate) * ( cos(theta) - cos(theta + noisy_yaw_rate * delta_t) );
			theta += noisy_yaw_rate * delta_t;
		}

		particles[i].x = x + d_x(gen);
		particles[i].y = y + d_y(gen);
		// particles[i].theta = normalize_angle(theta + d_theta(gen));
		particles[i].theta = theta + d_theta(gen);

	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

		// pre calculations for weight computation
		double gauss_norm = (1  / (2 * M_PI * std_landmark[0] * std_landmark[1]));
		double x_under = 2 * std_landmark[0] * std_landmark[0];
		double y_under = 2 * std_landmark[1] * std_landmark[1];

		// will be used for probability normalization
		double particle_weight_sum = 0.0;

		weights.clear();

		for (int i=0; i<num_particles; i++){  // loop over all particles
				Particle p = particles[i];
				p.weight = 1.0;
				p.sense_x.clear();
				p.sense_y.clear();
				p.associations.clear();

				for (int j=0; j<observations.size(); j++){  // loop over all observations
						LandmarkObs o = observations[j];

						// transform observation coordinates
						double x_m = p.x + cos(p.theta) * o.x - sin(p.theta) * o.y,
						       y_m = p.y + sin(p.theta) * o.x + cos(p.theta) * o.y;

						// save transformed observation in particle
						p.sense_x.push_back(x_m);
						p.sense_y.push_back(y_m);

						// look for closest landmark
						Map::single_landmark_s closest_landmark = map_landmarks.landmark_list[0];
						double min_dist = dist(x_m, y_m, closest_landmark.x_f, closest_landmark.y_f);
						Map::single_landmark_s l;
						double temp_dist;
						for (int k=0; k<map_landmarks.landmark_list.size(); k++){
								l = map_landmarks.landmark_list[k];
								temp_dist = dist(x_m, y_m, l.x_f, l.y_f);
								if (temp_dist < min_dist){
									closest_landmark = l;
									min_dist = temp_dist;
								}
						} // end for loop landmarks

						p.associations.push_back(closest_landmark.id_i);  // save id of closest landmark

						// compute weight of this observation
						double x_term = (x_m - closest_landmark.x_f) * (x_m - closest_landmark.x_f) / x_under;
						double y_term = (y_m - closest_landmark.y_f) * (y_m - closest_landmark.y_f) / y_under;
						double exponential = -(x_term + y_term);
						double w = gauss_norm * exp(exponential);

						double w_obsi = 1/(2*M_PI*std_landmark[0]*std_landmark[1])*exp(-0.5*
				    ((x_m - closest_landmark.x_f)*(x_m - closest_landmark.x_f)/(std_landmark[0]*std_landmark[0]) +
				    (y_m - closest_landmark.y_f)*(y_m - closest_landmark.y_f)/(std_landmark[1]*std_landmark[1])));

						// update total particle weight
						p.weight *= w_obsi;
				}  // end for loop observations

				// save particle weight and increment total particle weight sum
				weights.push_back(p.weight);
				particle_weight_sum += p.weight;

				// save updated particle
				particles[i] = p;
		} // end for loop particles

		// for(int i=0; i<weights.size(); i++){
		// 	weights[i] = weights[i] / particle_weight_sum;
		// }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::vector<Particle> new_particles;
	default_random_engine gen;

	// double particle_weight_sum = 0.0;
	// for (int i=0; i<particles.size(); i++)
	// 	particle_weight_sum += particles[i].weight;

	std::discrete_distribution<int> d(weights.begin(), weights.end());
	for(int i=0; i<particles.size(); i++){
		int sampled_particle_idx = d(gen);
		new_particles.push_back(particles[sampled_particle_idx]);
	}
	particles.clear();
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations = associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
