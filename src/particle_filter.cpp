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
		weights.push(1.0);
		Particle p;
		p.id = i;
		p.x = d_x(gen);
		p.y = d_y(gen);
		p.theta = d_theta(gen);
		p.weight = 1.0;
		particles.push(p);
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// std_pos = velocity and yaw rate measurement uncertainty
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	std::normal_distribution<double> d_velocity(velocity, std_pos[0]);
	std::normal_distribution<double> d_yaw_rate(yaw_rate, std_pos[1]);
	for(int i=0; i<num_particles; i++){
		double noisy_v = velocity + d_velocity(gen);
		double noisy_yaw_rate = yaw_rate + d_yaw_rate(gen);
		double new_theta = particles[i].theta + noisy_yaw_rate * delta_t;
		particles[i].theta = normalize_angle(new_theta);
		particles[i].x += noisy_v * cos(new_theta) * delta_t;
		particles[i].y += noisy_v * sin(new_theta) * delta_t;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	double temp_dist;
	for (int i=0; i<observations.size(); i++){
		int max_id = predicted[0].id;
		double max_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
		for (int j=1; j<predicted.size(); j++){
			temp_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (temp_dist < max_dist){
				max_id = predicted[j].id;
				max_dist = temp_dist;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// pre calculations for weight calculations
	double static first_term = (1  / (2 * M_PI * std_landmark[0] * std_landmark[1]));
	double static x_under = 2 * std_landmark[0] * std_landmark[0];
	double static y_under = 2 * std_landmark[0] * std_landmark[0];

	for (i=0; i<particles.size(); i++){
		// transform observations to map coordinates from particle perspective
		std::vector<LandmarkObs> transformedObservations;
		std::vector<LandmarkObs> landmarksInRange;
		for (j=0; j<observations.size(); j++){
			LandmarkObs tObs;
			tObs.x = particles[i].x +
							 cos(particles[i].theta) * observations[j].x -
							 sin(particles[i].theta) * observations[j].y;
		  tObs.y = particles[i].y +
							 sin(particles[i].theta) * observations[j].x +
							 cos(particles[i].theta) * observations[j].y;
			transformedObservations.push(tObs);
		} // end for transform observations

		// select map landmarks within the sensor range
		for (j=0; j<map_landmarks.landmark_list.size(); j++){
			if(dist(particle[i].x, article[i].y,
				 map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range){
				LandmarkObs tObs;
				tObs.x = map_landmarks.landmark_list[j].x_f;
				tObs.y = map_landmarks.landmark_list[j].y_f;
				// tObs.id = map_landmarks.landmark_list[j].id_i;
				tObs.id = j;
				landmarksInRange.push(tObs);
			} // end if
		}  // end for landmarks within range

		// associate each transformed measurement with a map landmark using NN
		dataAssociation(landmarksInRange, transformedObservations)

		// compute total weight my multiplying all weights of predicted landmaks
		// compute weight of each predicted landmark using a multivariate distribution
		double total_weight = 1;
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;
		for (j=0; j<transformedObservations.size(); j++){
			int associated_map_idx = transformedObservations[i].id;
			static double x = transformedObservations[i].x,
									  y = transformedObservations[i].y;
			static double mu_x = map_landmarks.landmark_list[associated_map_idx].x,
								    mu_y = map_landmarks.landmark_list[associated_map_idx].y;
			static double x_term = (x - mu_x) * (x - mu_x) / x_under;
			static double y_term = (y - mu_y) * (y - mu_y) / y_under;
			double w = first_term * exp(-(x_term + y_term));
			total_weight *= w;

			associations.push(map_landmarks.landmark_list[associated_map_idx].id_i);
			sense_x.push(x);
			sense_y.push(y);
		} // end for weights
		particles[i].weight = total_weight;

		SetAssociations(particles[i], associations, sense_x, sense_y);
		
	} //end for loop over particles

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
