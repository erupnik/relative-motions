#include "all.h"
#include "options.h"

namespace po = boost::program_options;


extern int cov_in_motions_main(InputFiles&,
                        LocalBundleOptions& ,
                        GlobalBundleOptions& ,
                        bool);


int main(int argc,char** argv)
{
    try
    {
        FLAGS_log_dir = "./";
        google::InitGoogleLogging(argv[0]);

        //input txt files
        InputFiles inputs;
        
        // bundle adjustment params
        LocalBundleOptions lba_opts ;
        GlobalBundleOptions gba_opts;

        // logging 
        int FLAGS_v;

        // internal checks 
        bool ceres_covariance=false;

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("tracks", po::value<std::string>(), "input files: feature tracks")
            ("views", po::value<std::string>(), "input files: local views")
            ("similarities", po::value<std::string>(), "input files: local2global frame similarities")
            ("global_poses_in", po::value<std::string>(), "input files: initial global poses")
            ("global_poses_out", po::value<std::string>(), "input files: output, refined global poses")

            ("lba_run_prop", po::value<bool>(), "local BA: covariance propagation ON/OFF")
            ("lba_feat_pds", po::value<double>(), "local BA: norm feature PDS")
            ("lba_R_pds", po::value<double>(), "local BA: rotation PDS")
            ("lba_C_pds", po::value<double>(), "local BA: center PDS")
            ("lba_loss_fea", po::value<double>(), "local BA: Huber loss threshold for features")
            ("lba_write_cov", po::value<bool>(), "local BA: covariance propagation ON/OFF")
            ("gba_R_pds", po::value<double>(), "global BA: rotation PDS")
            ("gba_C_pds", po::value<double>(), "global BA: center PDS")
            ("gba_loss_sim", po::value<double>(), "global BA: Huber loss threshold for similarities")
            ("gba_loss_gp", po::value<double>(), "global BA: Huber loss threshold for global poses")
            ("gba_inner_iter", po::value<bool>(), "global BA: Use inner iterations")
            
            ("ceres_cov", po::value<bool>(), "internal: compute ceres covariance")
;

       
        po::variables_map vmap;        
        po::store(po::parse_command_line(argc, argv, desc), vmap);
        po::notify(vmap);

        if (vmap.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }

        // input files 
        if (vmap.count("tracks"))
            inputs.tracks_file = vmap["tracks"].as<std::string>();
        else
        {
            throw "tracks not defined!";
            return 1.0;
        }
        if (vmap.count("views"))
            inputs.views_file = vmap["views"].as<std::string>();
        else
        {
            throw "views not defined!";
            return 1.0;
        }

        if (vmap.count("similarities"))
            inputs.similarities_file = vmap["similarities"].as<std::string>();
        else
        {
            throw "similarities not defined!";
            return 1.0;
        }
        if (vmap.count("global_poses_in"))
            inputs.global_poses_file = vmap["global_poses_in"].as<std::string>();
        else
        {
            throw "global poses not defined!";
            return 1.0;
        }
        if (vmap.count("global_poses_out"))
            inputs.output_poses_file = vmap["global_poses_out"].as<std::string>();
    
        if (vmap.count("ceres_cov"))
        {
            ceres_covariance = vmap["ceres_cov"].as<bool>();
        }

        //local bundle adjustment options 
        if (vmap.count("lba_run_prop")) 
        {
            lba_opts._RUN_PROP = vmap["lba_run_prop"].as<bool>() ;
        }
        if (vmap.count("lba_feat_pds")) 
        {
            lba_opts._FEAT_PDS = vmap["lba_feat_pds"].as<double>() ;
            //std::cout << "lba_feat_pds was set to " << lba_opts._FEAT_PDS <<  "\n";
        }
        
        if (vmap.count("lba_R_pds")) 
            lba_opts._ROT_PDS = vmap["lba_R_pds"].as<double>() ;

        if (vmap.count("lba_C_pds")) 
            lba_opts._C_PDS = vmap["lba_C_pds"].as<double>() ;

        if (vmap.count("lba_loss_fea")) 
            lba_opts._HUBER = vmap["lba_loss_fea"].as<double>() ;

        if (vmap.count("lba_write_cov"))
        {
            lba_opts._WRITE_COV = vmap["lba_write_cov"].as<bool>();
   
            // generate random key 
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> distrib(0, 100000);
            int random = distrib(gen);

            lba_opts._KEY += std::to_string(random);

        }

        //global bundle adjustemnt options
        if (vmap.count("gba_R_pds")) 
            gba_opts._ROT_PDS = vmap["gba_R_pds"].as<double>() ;

        if (vmap.count("gba_C_pds")) 
            gba_opts._C_PDS = vmap["gba_C_pds"].as<double>() ;

        if (vmap.count("gba_loss_sim")) 
            gba_opts._HUBER_S = vmap["gba_loss_sim"].as<double>() ;

        if (vmap.count("gba_loss_gp")) 
            gba_opts._HUBER_P = vmap["gba_loss_gp"].as<double>() ;

        if (vmap.count("gba_inner_iter")) 
            gba_opts._INNER_ITER = vmap["gba_inner_iter"].as<bool>() ;



        // run main program 
        cov_in_motions_main(inputs,lba_opts,gba_opts,
                            ceres_covariance);
        
    }
    catch (std::exception& e)
    {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) 
    {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;
}

