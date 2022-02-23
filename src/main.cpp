#include "all.h"


DEFINE_string(views_file,
              "",
              "The file containing relative motions");
DEFINE_string(tracks_file,
              "",
              "The file containing tracks (image features)");
DEFINE_string(similarities_file,
              "",
              "The file containing the 3D similarities between the relative and global poses");
DEFINE_string(global_poses_file,
              "",
              "The file containing the inital global poses");
DEFINE_string(output_poses_file,
              "output_poses.txt",
              "Compute covariance matrices per motions (this is only a check)");
DEFINE_bool(covariance,
              false,
              "Compute covariance matrices per motions (this is only a check)");


extern int cov_in_motions_main(std::string,
                               std::string,
                               std::string,
                               std::string,
                               std::string,
                               bool);


int main(int argc, char** argv)
{

  FLAGS_log_dir = "./";

  google::InitGoogleLogging(argv[0]);
  GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  LOG(INFO) << "cov_in_motions --views_file=" << FLAGS_views_file
            << " --tracks_file=" << FLAGS_tracks_file
            << " --similarities_file=" << FLAGS_similarities_file
            << " --global_poses_file=" << FLAGS_global_poses_file
            << " --output_poses_file=" << FLAGS_output_poses_file
            << " --covariance=" << FLAGS_covariance;


  return cov_in_motions_main(FLAGS_views_file,
                             FLAGS_tracks_file,
                             FLAGS_similarities_file,
                             FLAGS_global_poses_file,
                             FLAGS_output_poses_file,
                             FLAGS_covariance);


}
