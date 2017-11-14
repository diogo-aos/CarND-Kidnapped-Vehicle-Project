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
	// std_pos = velocity and yaw rate measurement uncertainty
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	std::normal_distribution<double> d_velocity(velocity, std_pos[0]);
	std::normal_distribution<double> d_yaw_rate(yaw_rate, std_pos[1]);
	for(int i=0; i<num_particles; i++){
		double x = particles[i].x,
				   y = particles[i].y,
					 theta = particles[i].theta;

		double noisy_v = velocity + d_velocity(gen);
		double noisy_yaw_rate = yaw_rate + d_yaw_rate(gen);

		if (fabs(yaw_rate) < 0.0001){  // yaw rate is 0
			x += noisy_v * delta_t * cos(theta);
			y += noisy_v * delta_t * sin(theta);
			// theta = theta;
		}
		else {  // yaw rate different than zero
			x += (noisy_v / noisy_yaw_rate) * ( sin(theta + noisy_yaw_rate * delta_t) - sin(theta) );
			y += (noisy_v / noisy_yaw_rate) * ( cos(theta) - cos(theta + noisy_yaw_rate * delta_t) );
			theta += noisy_yaw_rate * delta_t;
		}

		particles[i].x = x;
		particles[i].y = y;
		particles[i].theta = normalize_angle(theta);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	double temp_dist;

	for (int i=0; i<observations.size(); i++){
		double x_obs = observations[i].x,
		       y_obs = observations[i].y;
#if DEBUG == 1
		cout << "        dataAssociation :: observation " << i << endl;
#endif
		int associated_id = predicted[0].id;
		double max_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);

		for (int j=1; j<predicted.size(); j++){
			double x_pred = predicted[j].x,
			       y_pred = predicted[j].y;
			temp_dist = dist(x_obs, y_obs, x_pred, y_pred);
			if (temp_dist < max_dist){  // select predicted landmark closest to observation
				associated_id = predicted[j].id;
				max_dist = temp_dist;
			}
		}
#if DEBUG == 1
		cout << "        dataAssociation :: associated to " << associated_id << endl;
#endif
		observations[i].id = associated_id;
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

	// pre calculations for weight computation
	double gauss_norm = (1  / (2 * M_PI * std_landmark[0] * std_landmark[1]));
	double x_under = 2 * std_landmark[0] * std_landmark[0];
	double y_under = 2 * std_landmark[1] * std_landmark[1];

	double particle_weight_sum = 0.0;  // will be used for probability normalization

	for (int i=0; i<num_particles; i++){  // loop over all particles
		// load particle x, y, theta to more convenient variables
		double x_p = particles[i].x,
				   y_p = particles[i].y,
					 theta_p = particles[i].theta;
#if DEBUG_2 == 1
		cout << "    update :: particle " << i << " x_p=" << x_p << "  y_p=" << y_p << endl;
#endif

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		//  								TRANSFORM OBSERVATIONS
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#if DEBUG_2 == 1
		cout << "    update :: transformed observations" << endl;
#endif
		std::vector<LandmarkObs> transformedObservations;
		for (int j=0; j<observations.size(); j++){
			double x_c = observations[j].x,
			       y_c = observations[j].y;

			double x_m = x_p + cos(theta_p) * x_c - sin(theta_p) * y_c,
			       y_m = y_p + sin(theta_p) * x_c + cos(theta_p) * y_c;

			LandmarkObs tObs;
			tObs.x = x_m;
		  tObs.y = y_m;
			transformedObservations.push_back(tObs);
		} // end for transform observations

		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		//  								SELECT PREDICTED LANDMARKS
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#if DEBUG_2 == 1
		cout << "    update :: landmarkd in sensor range" << endl;
#endif
		std::vector<LandmarkObs> landmarksInRange;
		for (int j=0; j<map_landmarks.landmark_list.size(); j++){
			double x_l = map_landmarks.landmark_list[j].x_f,
			       y_l = map_landmarks.landmark_list[j].y_f;
			// check if distance between particle and landmark is within sensor range
			if(dist(x_p, y_p, x_l, y_l) <= sensor_range){
				LandmarkObs tObs;
				tObs.x = x_l;
				tObs.y = y_l;
				// tObs.id = map_landmarks.landmark_list[j].id_i;
				tObs.id = j;
				landmarksInRange.push_back(tObs);
			} // end if dist below sensor range
		}  // end for landmarks within range
#if DEBUG_2 == 1
		cout << "    update :: num predicted landmarks in sensor range " << landmarksInRange.size() << endl;
		cout << "    update :: num observations " << transformedObservations.size() << endl;
#endif

		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		//  								ASSOCIATE OBSERVATIONS WITH LANDMARKS
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#if DEBUG_2 == 1
		cout << "    update :: dataAssociation" << endl;
#endif
		dataAssociation(landmarksInRange, transformedObservations);

		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		//  								COMPUTE WEIGHTS
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#if DEBUG_2 == 1
		cout << "    update :: weight computation" << endl;
#endif
		double total_weight = 1.0;
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;
		for (int j=0; j<transformedObservations.size(); j++){
			// index of landmark associated with this observation
			int associated_map_idx = transformedObservations[j].id;

		  double x = transformedObservations[j].x,
						        y = transformedObservations[j].y;
			double mu_x = map_landmarks.landmark_list[associated_map_idx].x_f,
						        mu_y = map_landmarks.landmark_list[associated_map_idx].y_f;

			double x_term = (x - mu_x) * (x - mu_x) / x_under;
			double y_term = (y - mu_y) * (y - mu_y) / y_under;
			double exponential = -(x_term + y_term);
			double w = gauss_norm * exp(exponential);
			total_weight *= w;

#if DEBUG == 1
			cout <<	std::scientific;
			cout << "        update :: particle " << i << " observation " << j << endl;
			cout << "        update :: x_p = " << x_p << "  y_p = " << y_p << endl;
			cout << "        update :: x = " << x << "  y = " << y << endl;
			cout << "        update :: mu_x = " << mu_x << "  mu_y = " << mu_y << endl;
			cout << "        update :: x - mu_x = " << x - mu_x << endl;
			cout << "        update :: y - mu_y = " << y - mu_y << endl;
			cout << "        update :: x_term = " << x_term << endl;
			cout << "        update :: x_under = " << x_under << endl;
			cout << "        update :: gauss_norm = " << gauss_norm << endl;
			cout << "        update :: exponential = " << exponential << endl;
			cout << "        update :: w = " << w << endl;
			cout << "        update :: total_weight = " << total_weight << endl;
#endif

			associations.push_back(map_landmarks.landmark_list[associated_map_idx].id_i);
			sense_x.push_back(x);
			sense_y.push_back(y);
		} // end for weights
		particles[i].weight = total_weight;
		weights[i] = total_weight;
		particle_weight_sum += total_weight;

		particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);

	} //end for loop over particles

	// normalize particle weights/probabilities
	double normalized_w;
	double total_prob = 0.0;
	for (int i=0; i<weights.size(); i++){
		normalized_w = weights[i] / particle_weight_sum;
		weights[i] = normalized_w;
		particles[i].weight = normalized_w;
		total_prob += weights[i];
	}  // end for weights
} // end function

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
