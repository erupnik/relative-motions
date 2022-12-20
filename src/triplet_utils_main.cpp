#include "triplet_utils.h"

DEFINE_string(views_file,
              "",
              "The file containing relative motions");
DEFINE_string(similarities_file,
              "",
              "The file containing the 3D similarities between the relative and global poses");
DEFINE_string(global_poses_file,
              "",
              "The file containing the inital global poses");
DEFINE_string(output_poses_file,
              "output_poses.txt",
              "Output file with inlier (inlier_out_poses.txt) or outlier poses (outlier_out_poses.txt) or all of it output_poses.txt");
DEFINE_int32(filter,
              false,
              "Filter relative motions");
DEFINE_bool(do_only_predict,
              false,
              "Do only prediction");
DEFINE_string(filtered_view_file,
              "EGs_filtered.txt",
              "The file with inlier and outlier relative motions");

int main(int argc,char** argv)
{

   FLAGS_log_dir = "./";

   google::InitGoogleLogging(argv[0]);
   GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

   LOG(INFO) << "triplet_utils --views_file=" << FLAGS_views_file
            << " --similarities_file=" << FLAGS_similarities_file
            << " --global_poses_file=" << FLAGS_global_poses_file
            << " --output_poses_file=" << FLAGS_output_poses_file
            << " --filter=" << FLAGS_filter
            << " --filtered_file=" << FLAGS_filtered_view_file
            << " --do_only_predict=" << FLAGS_do_only_predict;
     
//    cTripletSet TriSet(FLAGS_views_file,FLAGS_similarities_file,FLAGS_global_poses_file);

    std::cout << "reussite!" << std::endl;

    cFilterTrip filter(FLAGS_views_file,
                       FLAGS_similarities_file,
                       FLAGS_do_only_predict,
                       FLAGS_filter,
                       FLAGS_output_poses_file,
                       FLAGS_global_poses_file,
                       FLAGS_filtered_view_file);
    /*cAppCovInMotion anAp(FLAGS_views_file,"",
                         FLAGS_similarities_file,
                         FLAGS_global_poses_file,
                         FLAGS_output_poses_file,
                         false);
  */ 
    return EXIT_SUCCESS;
}
