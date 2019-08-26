/*
 * Copyright 2018 The Cartographer Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "cartographer/ground_truth/autogenerate_ground_truth.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include "cartographer/mapping/proto/trajectory.pb.h"
#include "cartographer/transform/transform.h"
#include "glog/logging.h"

namespace cartographer {
namespace ground_truth {
namespace {

struct Error {
  double translation;
  double rotation;
  int node_trajectory_id;
};

void MeanAndStdDev(const std::vector<double>& values, double& mean, double& standard_deviation) {
  mean =
      std::accumulate(values.begin(), values.end(), 0.) / values.size();
  double sum_of_squared_differences = 0.;
  for (const double value : values) {
    sum_of_squared_differences += common::Pow2(value - mean);
  }
  standard_deviation =
      std::sqrt(sum_of_squared_differences / (values.size() - 1));

}

void WriteStatisticsToFile(const std::vector<Error>& local_errors,
                          const std::vector<Error>& global_errors,
                          double number_of_trajectories,
                           std::ofstream& statistics_errors_file) {

  statistics_errors_file
      << "trajectory_id, local_loop_closure_count, local_translational_mean_error, local_translational_std_error,"
         "local_rotational_std_error, local_rotational_std_error,"
         "global_loop_closure_count, global_translational_mean_error, global_translational_std_error,"
         "global_rotational_std_error, global_rotational_std_error\n";

  std::vector<std::vector<double>> local_translational_errors(number_of_trajectories);
  std::vector<std::vector<double>> local_rotational_errors_degrees(number_of_trajectories);

  std::vector<std::vector<double>> global_translational_errors(number_of_trajectories);
  std::vector<std::vector<double>> global_rotational_errors_degrees(number_of_trajectories);

  for (const Error& error : local_errors) {
    local_translational_errors[error.node_trajectory_id].push_back(error.translation);
    local_rotational_errors_degrees[error.node_trajectory_id].push_back(
        common::RadToDeg(error.rotation));
  }

  for (const Error& error : global_errors) {
    global_translational_errors[error.node_trajectory_id].push_back(error.translation);
    global_rotational_errors_degrees[error.node_trajectory_id].push_back(
        common::RadToDeg(error.rotation));
  }
  
  for(int id = 0; id < number_of_trajectories; ++id) {
    statistics_errors_file << id << ",";
    statistics_errors_file << local_translational_errors[id].size() << ",";
    double local_trans_mean = 0.0, local_trans_std = 0.0;
    MeanAndStdDev(local_translational_errors[id], local_trans_mean, local_trans_std);
    statistics_errors_file <<  local_trans_mean << ","
                         << local_trans_std << ",";
           
    double local_rotational_mean = 0.0, local_rotational_std = 0.0;      
    MeanAndStdDev(local_rotational_errors_degrees[id], local_rotational_mean, local_rotational_std);
    statistics_errors_file <<  local_rotational_mean << ","
                           << local_rotational_std << ",";
    
    statistics_errors_file << global_translational_errors[id].size() << ",";
    double global_trans_mean = 0.0, global_trans_std = 0.0;
    MeanAndStdDev(global_translational_errors[id], global_trans_mean, global_trans_std);
    statistics_errors_file <<  global_trans_mean << ","
                         << global_trans_std << ",";
           
    double global_rotational_mean = 0.0, global_rotational_std = 0.0;      
    MeanAndStdDev(global_rotational_errors_degrees[id], global_rotational_mean, global_rotational_std);
    statistics_errors_file <<  global_rotational_mean << ","
                           << global_rotational_std << "\n";
  }
         
}

void WriteRelationMetricsToFile(const Error& error,
                                const transform::Rigid3d& pose,
                                const int64& timestamp,
                                const double covered_distance,
                                const int node_trajectory_id,
                                const int submap_trajectory_id,
                                std::ofstream& relation_errors_file
                                ) {

    long int epoch = 621355968000000000;
    double time = (timestamp- epoch)/pow(10,7);
    double translational_error = error.translation;
    double rotational_errors_degree =
        common::RadToDeg((error.rotation));
    relation_errors_file << translational_error << ","
                         << rotational_errors_degree << ","
                         << pose.translation().x() << ","
                         << pose.translation().y() << ","
                         << transform::GetYaw(pose) << ","
                         << std::setprecision (std::numeric_limits<double>::digits10 + 1) << time << ","
                         << covered_distance << ","
                         << node_trajectory_id << ","
                         << submap_trajectory_id << "\n";
  
}

std::vector<double> ComputeCoveredDistance(
    const mapping::proto::Trajectory& trajectory) {
  std::vector<double> covered_distance;
  covered_distance.push_back(0.);
  CHECK_GT(trajectory.node_size(), 0)
      << "Trajectory does not contain any nodes.";
  for (int i = 1; i < trajectory.node_size(); ++i) {
    const auto last_pose = transform::ToRigid3(trajectory.node(i - 1).pose());
    const auto this_pose = transform::ToRigid3(trajectory.node(i).pose());
    covered_distance.push_back(
        covered_distance.back() +
        (last_pose.inverse() * this_pose).translation().norm());
  }
  return covered_distance;
}

// We pick the representative node in the middle of the submap.
//
// TODO(whess): Should we consider all nodes inserted into the submap and
// exclude, e.g. based on large relative linear or angular distance?
std::vector<int> ComputeSubmapRepresentativeNode(
    const mapping::proto::PoseGraph& pose_graph, int trajectory_id) {
  std::vector<int> submap_to_node_index;
  for (const auto& constraint : pose_graph.constraint()) {
    if (constraint.tag() != mapping::proto::PoseGraph::Constraint::INTRA_SUBMAP || 
        constraint.submap_id().trajectory_id() != trajectory_id ||
        constraint.node_id().trajectory_id() != trajectory_id) {
      continue;
    } 

    CHECK_EQ(constraint.submap_id().trajectory_id(), trajectory_id);
    CHECK_EQ(constraint.node_id().trajectory_id(), trajectory_id);

    const int next_submap_index = static_cast<int>(submap_to_node_index.size());
    const int submap_index = constraint.submap_id().submap_index();
    
    if (submap_index <= next_submap_index) {
      continue;
    }

    CHECK_EQ(submap_index, next_submap_index + 1);
    submap_to_node_index.push_back(constraint.node_id().node_index());
  }
  return submap_to_node_index;
}

}  // namespace

proto::GroundTruth GenerateGroundTruthForMultipleTrajectories(
    const mapping::proto::PoseGraph& pose_graph,
    const double min_covered_distance, const double outlier_threshold_meters,
    const double outlier_threshold_radians, const std::string& output_filename) {
  
  const std::string relation_metrics_filename = output_filename + "_relations.csv";
  const std::string statistics_error_filename = output_filename + "_statistics.csv";

  proto::GroundTruth ground_truth;
  std::ofstream relation_errors_file;
  relation_errors_file.open(relation_metrics_filename);
  relation_errors_file << "translational_error,rotational_errors_degree,"
      "pose_x,pose_y,pose_yaw,timestamp,covered_distance,node_trajectory_id,submap_trajectory_id\n";

  std::ofstream stats_errors_file;
  stats_errors_file.open(statistics_error_filename);
  
  LOG(INFO) << "Writing relation metrics to '" + relation_metrics_filename +
                   "'...";
  

  std::vector<std::vector<int> >trajectory_submap_node_index;
  std::vector<std::vector<double> >trajectory_covered_distance;
  std::vector<Error> local_error_vec;
  std::vector<Error> global_error_vec;
  int number_of_trajectories = pose_graph.trajectory().size();
  
  // compute submap representative node for all trajectories
  for (size_t trajectory_id = 0; trajectory_id < number_of_trajectories; ++trajectory_id) {
    const mapping::proto::Trajectory& trajectory = pose_graph.trajectory(trajectory_id);
    const std::vector<double> covered_distance =
        ComputeCoveredDistance(trajectory);
    trajectory_covered_distance.push_back(covered_distance);

    const std::vector<int> submap_to_node_index =
        ComputeSubmapRepresentativeNode(pose_graph, trajectory_id);
    trajectory_submap_node_index.push_back(submap_to_node_index);
  }

  int num_outliers = 0;
  int num_relations = 0;
  for (const auto& constraint : pose_graph.constraint()) {
    // We're only interested in local loop closure constraints.
    if (constraint.tag() ==
        mapping::proto::PoseGraph::Constraint::INTRA_SUBMAP ) {
      continue;
    }
    int node_trajectory_id = constraint.node_id().trajectory_id();
    int submap_trajectory_id = constraint.submap_id().trajectory_id();
    const mapping::proto::Trajectory& node_trajectory = pose_graph.trajectory(node_trajectory_id);
    const mapping::proto::Trajectory& submap_trajectory = pose_graph.trajectory(submap_trajectory_id);
    const std::vector<int>& submap_to_node_index = trajectory_submap_node_index[submap_trajectory_id];
    const std::vector<double>& covered_distance_trajectory_node = trajectory_covered_distance[node_trajectory_id];
    const std::vector<double>& covered_distance_trajectory_submap = trajectory_covered_distance[submap_trajectory_id];

    if (constraint.submap_id().submap_index() >=
          static_cast<int>(submap_to_node_index.size())) {
        continue;
    }
    
    const int matched_node = constraint.node_id().node_index();
    const int representative_node =
          submap_to_node_index.at(constraint.submap_id().submap_index());
    // Compute the transform between the nodes according to the solution and
    // the constraint.
    const transform::Rigid3d solution_pose1 =
        transform::ToRigid3(submap_trajectory.node(representative_node).pose());
    const transform::Rigid3d solution_pose2 =
        transform::ToRigid3(node_trajectory.node(matched_node).pose());
    const transform::Rigid3d solution =
        solution_pose1.inverse() * solution_pose2;

    const transform::Rigid3d submap_solution = transform::ToRigid3(
        submap_trajectory.submap(constraint.submap_id().submap_index()).pose());
    const transform::Rigid3d submap_solution_to_node_solution =
        solution_pose1.inverse() * submap_solution;
    const transform::Rigid3d node_to_submap_constraint =
        transform::ToRigid3(constraint.relative_pose());
    const transform::Rigid3d expected =
        submap_solution_to_node_solution * node_to_submap_constraint;
    
    const transform::Rigid3d error = solution * expected.inverse();
    Error error_struct({error.translation().norm(),
             (transform::GetAngle(error)), node_trajectory_id});

    double covered_distance_in_constraint = 
        std::abs(covered_distance_trajectory_node.at(matched_node) -
                 covered_distance_trajectory_submap.at(representative_node));
    
    int start_trajectory_id = submap_trajectory_id > node_trajectory_id ? node_trajectory_id : submap_trajectory_id;
    int end_trajectory_id = submap_trajectory_id > node_trajectory_id ? submap_trajectory_id : node_trajectory_id;
    
    for (int trajectory_id = (start_trajectory_id); trajectory_id < end_trajectory_id-1; ++trajectory_id) {
      covered_distance_in_constraint += trajectory_covered_distance[trajectory_id].back();
    }
    if (covered_distance_in_constraint < min_covered_distance) {
      continue;
    }

    if(node_trajectory_id == submap_trajectory_id) {
      local_error_vec.push_back(error_struct);
    } else if(node_trajectory_id > submap_trajectory_id){
      global_error_vec.push_back(error_struct);
    }

    WriteRelationMetricsToFile(error_struct, 
                          solution_pose2, 
                          node_trajectory.node(matched_node).timestamp(), 
                          covered_distance_in_constraint, 
                          node_trajectory_id, 
                          submap_trajectory_id, 
                          relation_errors_file);
    
    if (error.translation().norm() > outlier_threshold_meters ||
        transform::GetAngle(error) > outlier_threshold_radians) {;
      ++num_outliers;
      continue;
    }

    auto* const new_relation = ground_truth.add_relation();
    new_relation->set_timestamp1(
        submap_trajectory.node(representative_node).timestamp());
    new_relation->set_timestamp2(node_trajectory.node(matched_node).timestamp());
    *new_relation->mutable_expected() = transform::ToProto(expected);
    new_relation->set_covered_distance(covered_distance_in_constraint);
    new_relation->set_trajectory_id(constraint.submap_id().trajectory_id());
  }
  WriteStatisticsToFile(local_error_vec, global_error_vec, number_of_trajectories, stats_errors_file);
  relation_errors_file.close();  
  stats_errors_file.close();
  return ground_truth;
}

proto::GroundTruth GenerateGroundTruth(
    const mapping::proto::PoseGraph& pose_graph,
    const double min_covered_distance, const double outlier_threshold_meters,
    const double outlier_threshold_radians) {
  const mapping::proto::Trajectory& trajectory = pose_graph.trajectory(0);
  const std::vector<double> covered_distance =
      ComputeCoveredDistance(trajectory);

  const std::vector<int> submap_to_node_index =
      ComputeSubmapRepresentativeNode(pose_graph, 0);

  int num_outliers = 0;
  proto::GroundTruth ground_truth;
  for (const auto& constraint : pose_graph.constraint()) {
    // We're only interested in loop closure constraints.
    if (constraint.tag() ==
        mapping::proto::PoseGraph::Constraint::INTRA_SUBMAP) {
      continue;
    }

    // For some submaps at the very end, we have not chosen a representative
    // node, but those should not be part of loop closure anyway.
    CHECK_EQ(constraint.submap_id().trajectory_id(), 0);
    CHECK_EQ(constraint.node_id().trajectory_id(), 0);
    if (constraint.submap_id().submap_index() >=
        static_cast<int>(submap_to_node_index.size())) {
      continue;
    }
    const int matched_node = constraint.node_id().node_index();
    const int representative_node =
        submap_to_node_index.at(constraint.submap_id().submap_index());

    // Covered distance between the two should not be too small.
    double covered_distance_in_constraint =
        std::abs(covered_distance.at(matched_node) -
                 covered_distance.at(representative_node));
    if (covered_distance_in_constraint < min_covered_distance) {
      continue;
    }

    // Compute the transform between the nodes according to the solution and
    // the constraint.
    const transform::Rigid3d solution_pose1 =
        transform::ToRigid3(trajectory.node(representative_node).pose());
    const transform::Rigid3d solution_pose2 =
        transform::ToRigid3(trajectory.node(matched_node).pose());
    const transform::Rigid3d solution =
        solution_pose1.inverse() * solution_pose2;

    const transform::Rigid3d submap_solution = transform::ToRigid3(
        trajectory.submap(constraint.submap_id().submap_index()).pose());
    const transform::Rigid3d submap_solution_to_node_solution =
        solution_pose1.inverse() * submap_solution;
    const transform::Rigid3d node_to_submap_constraint =
        transform::ToRigid3(constraint.relative_pose());
    const transform::Rigid3d expected =
        submap_solution_to_node_solution * node_to_submap_constraint;

    const transform::Rigid3d error = solution * expected.inverse();

    if (error.translation().norm() > outlier_threshold_meters ||
        transform::GetAngle(error) > outlier_threshold_radians) {
      ++num_outliers;
      continue;
    }
    auto* const new_relation = ground_truth.add_relation();
    new_relation->set_timestamp1(
        trajectory.node(representative_node).timestamp());
    new_relation->set_timestamp2(trajectory.node(matched_node).timestamp());
    *new_relation->mutable_expected() = transform::ToProto(expected);
    new_relation->set_covered_distance(covered_distance_in_constraint);
  }
  LOG(INFO) << "Generated " << ground_truth.relation_size()
            << " relations and ignored " << num_outliers << " outliers.";
  return ground_truth;
}


}  // namespace ground_truth
}  // namespace cartographer
