#include "all.h"
#include "cov_in_motions.h"

// TODO
//   -
//
// to consider :
//     - sparse matrices for H

/*

  3 relative poses: i, i in [1..3]
     - H_i = (J_t W J_)    6 x 6  hessian
     - g_i = J_t W res     6 x 1  gradient vector

  1 similarity: s, alpha, beta
       pose_i = s * alpha * Pose_i + beta //global to local


  Jacobians, hessians:
    Hi =  J_rt_RT^t @  H_i  @ J_rt_RT
    gi =  J_rt_RT^t @ g_i

    H = SUM |         |   => can become big for Nk images so think of sparse matrices
            |   Hi    |
            |         |

    where,
    J_rt_RT = [ J_rt_R | 0      ]    jacobian
              [   0    | J_rt_T ]
    J_rt_R =  [ alpha @ [1,1,1]x ]
           =   |  a11 a12  a13 |
               |  a21 a22  a23 | @ [1,1,1]x
               |  a31 a32  a33 |
    J_rt_T =  [ s * alpha  ]
                  |  a11 a12  a13 |
          =    s *|  a21 a22  a23 |
                  |  a31 a32  a33 |

    sizes for M=3 images:
      Hi        6 x  6
      gi        6 x  1
      H        18 x 18
      g        18 x  1
      J_rt_RT

  Solution:
     delta =  H^-1 @ g    18 x 1

  Update: Pose+ = exp(delta)^T * Pose
        [R+,T+] = exp(delta)^T * [R | T]

*/
Mat3d skew_symmetric(Vec3d w)
{
  Mat3d ss_mat;
  ss_mat << 0, -w[2], w[1],
            w[2], 0, -w[0],
           -w[1], w[0], 0;

  return ss_mat;
}

void test_gen_2motions(cNviewPoseT<Mat6d,Vec6d>*& motion_1,cNviewPoseT<Mat12d,Vec12d>*& motion_2)
{
    //relative poses
    //motion1
    Mat3d r3_m1 = Mat3d::Identity();
    Vec3d t3_m1 = Vec3d::Zero();

    Mat3d r4_m1;
    r4_m1 <<  0.990304, -0.0184962, -0.137682,
              0.0236068, 0.999088, 0.0355792,
              0.136898, -0.0384845, 0.989837;
    Vec3d t4_m1;
    t4_m1 <<  0.998794, -0.0280127, 0.0403243;

    //relative poses
    cPoseBasic*  pose3_m1 = new cPoseBasic(r3_m1,t3_m1,"0003.png");
    cPoseBasic*  pose4_m1 = new cPoseBasic(r4_m1,t4_m1,"0004.png");

    //hessian and gradient
    Mat6d hessian2;
    hessian2 << 0, 0, 0, 0, 0, 0,
               0, 136.096, 0.655852, 274.055, -11.6732, 23.0635,
               0, 0.655852, 106.606, 4.26845, -69.0993, -2.56514,
               0, 274.055, 4.26845, 88743.8, -50.4043, 103.143,
               0, -11.6732,-69.0993, -50.4043, 89113.6, 99.535,
               0, 23.0635, -2.56514, 103.143, 99.535, 86756.9;
    Vec6d gradient2;
    gradient2 << 0, 0.0829167, 0.021951, 0.738334, -0.0964338, 0.18924;//0.18924;

    //cHessianGradient2 * aHg2 = new cHessianGradient2(hessian2,gradient2);
    cHessianGradientN<Mat6d,Vec6d> * aHgN = new cHessianGradientN<Mat6d,Vec6d>(hessian2,gradient2);
    //aHg->printH();

    //affine transformation
    Mat3d alpha_m1;
    alpha_m1 <<  -0.943879, 0.0336532, 0.328572,
                 0.0452759, 0.998588, 0.0277849,
                 -0.327173, 0.041102, -0.94407;
    Mat3d alpha_m1_inv = alpha_m1.inverse();
    Vec3d beta_m1;
    beta_m1 << 1.45337, 0.782566, 4.01418;
    double s_m1 = 2.50449;

    motion_1 = (new cNviewPoseT<Mat6d,Vec6d>(pose3_m1,pose4_m1,NULL,aHgN));
    motion_1->lambda() = s_m1;
    motion_1->alpha() = alpha_m1;
    motion_1->beta() = beta_m1;

    motion_1->PrintAlpha();
    motion_1->Hg_().printH();

    //motion2
    Mat3d r2_m2 = Mat3d::Identity();
    Vec3d t2_m2 = Vec3d::Zero();

    Mat3d r3_m2;
    r3_m2 <<  0.968305, -0.0471622, -0.245278,
              0.0500244, 0.998733, 0.00544876,
              0.24471, -0.0175459, 0.969438;
    Vec3d t3_m2;
    t3_m2 <<  0.294887, 0.0851058, 0.36799;

    Mat3d r4_m2;
    r4_m2 << 0.923964, -0.055705, -0.378401,
             0.073779, 0.996714, 0.0334227,
             0.375296, -0.0587994, 0.925038;
    Vec3d t4_m2;
    t4_m2 <<  0.841244, 0.100972, 0.531143;

    //relative poses
    cPoseBasic*  pose2_m2 = new cPoseBasic(r2_m2,t2_m2,"0002.png");
    cPoseBasic*  pose3_m2 = new cPoseBasic(r3_m2,t3_m2,"0003.png");
    cPoseBasic*  pose4_m2 = new cPoseBasic(r4_m2,t4_m2,"0004.png");

    //hessian and gradient
    Mat12d hessian3;
    hessian3 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 157.704, 0.94878, 245.912, -15.638, -0.599963, 0, 0, 0, 0, 0, 0,
               0, 0.94878, 108.609, 1.595, -61.9846, -1.17501, 0, 0, 0, 0, 0, 0,
               0, 245.912, 1.595, 87674.9, -22.8704, -26.9289, 0, 0, 0, 0, 0, 0,
               0,-15.638, -61.9846, -22.8704, 87786.9, 22.552, 0, 0, 0, 0, 0, 0,
               0, -0.599963, -1.17501, -26.9289, 22.552, 86590.7, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 148.867, -0.395078, 18.2715, -16.9635, -242.946, -10.8299,
               0, 0, 0, 0, 0, 0, -0.395078, 156.773, 1.03086, 245.458, -20.8779, -7.17429,
               0, 0, 0, 0, 0, 0, 18.2715, 1.03086, 113.777, -1.24273, -90.2474, -4.1709,
               0, 0, 0, 0, 0, 0, -16.9635, 245.458, -1.24273, 87679.6, -13.2412, -57.8773,
               0, 0, 0, 0, 0, 0, -242.946, -20.8779, -90.2474, -13.2412, 87810.4, 59.2815,
               0, 0, 0, 0, 0, 0, -10.8299, -7.17429, -4.1709, -57.8773, 59.2815, 86603;
    Vec12d gradient3;
    gradient3 << 0, -0.0593826, -0.297001,
                -0.0269822, 0.569161, -0.0778623,
                0.227898, 0.103299, 0.344239,
                0.304507, -0.314465, -0.0427875;


    cHessianGradientN<Mat12d,Vec12d> * aHgN_ =
                     new cHessianGradientN<Mat12d,Vec12d>(hessian3,gradient3);

    //affine transformation
    Mat3d alpha_m2;
    alpha_m2 <<  -0.835362, -0.0244854, 0.549154,
                 -0.00363671, 0.999232, 0.0390211,
                -0.549688, 0.0305996, -0.834809;
    Mat3d alpha_m2_inv = alpha_m2.inverse();
    Vec3d beta_m2;
    beta_m2 << 0.510063, 0.616558, 2.93794;
    double s_m2 = 1.52224;


    motion_2 = (new cNviewPoseT<Mat12d,Vec12d>(pose2_m2,pose3_m2,pose4_m2,aHgN_));
    motion_2->lambda() = s_m2;
    motion_2->alpha() = alpha_m2;
    motion_2->beta() = beta_m2;

    motion_2->PrintAlpha();
    motion_2->Hg_().printH();

    std::cout << "above H" ;
    motion_1->View(1).Show();
    motion_2->View(2).Show();

}

void test_propagate_cov(cNviewPoseT<Mat6d,Vec6d>*& motion_1)
{
  Mat3d J_rt_T = motion_1->lambda() * motion_1->alpha();
  Mat3d skew_sym = skew_symmetric(Vec3d::Ones());
  Mat3d J_rt_R = motion_1->alpha() * skew_sym;

  Mat6d J_rt_RT = Mat6d::Zero();
  J_rt_RT.block(0,0,3,3) = J_rt_T;
  J_rt_RT.block(3,3,3,3) = J_rt_R;

  //std::cout << "\nJ_rt_RT\n" << J_rt_RT << '\n';
  //std::cout << "\nmotion_1->Hg_().H_()\n" << motion_1->Hg_().H_() << '\n';

  Mat6d JtHJ = J_rt_RT.transpose() * motion_1->Hg_().H_() * J_rt_RT;
  //std::cout << "\nJtHJ\n" << JtHJ << '\n';
  motion_1->Hg_().H_() = JtHJ;
  motion_1->Hg_().g_() = J_rt_RT.transpose() * motion_1->Hg_().g_();
  std::cout << "\nmotion_1->Hg_().H_()\n" << motion_1->Hg_().H_() << '\n';
  std::cout << "\n =============================\n";
}

void test_propagate_cov(cNviewPoseT<Mat12d,Vec12d>*& motion)
{
  Mat3d skew_sym = skew_symmetric(Vec3d::Ones());
  Mat3d J_rt_T = motion->lambda() * motion->alpha();
  Mat3d J_rt_R = motion->alpha() * skew_sym;

  Mat12d J_rt_RT = Mat12d::Zero();
  J_rt_RT.block(0,0,3,3) = J_rt_T;
  J_rt_RT.block(3,3,3,3) = J_rt_R;
  J_rt_RT.block(6,6,3,3) = J_rt_T;
  J_rt_RT.block(9,9,3,3) = J_rt_R;

  //std::cout << "\nJ_rt_RT\n" << J_rt_RT << '\n';
  //std::cout << "\nmotion->Hg_().H_()\n" << motion->Hg_().H_() << '\n';

  Mat12d JtHJ = J_rt_RT.transpose() * motion->Hg_().H_() * J_rt_RT;
  //std::cout << "\nJtHJ\n" << JtHJ << '\n';
  motion->Hg_().H_() = JtHJ;
  motion->Hg_().g_() = J_rt_RT.transpose() * motion->Hg_().g_();
  std::cout << "\nmotion->Hg_().H_()\n" << motion->Hg_().H_() << '\n';
  std::cout << "\n =============================\n";
}

void test_similarity_trafo()
{
    //global to local transformation :
    //pose = s * alpha * Pose + beta

    //local to global transformation :
    // Pose = 1/s * alpha^-1 * (pose - beta)

    /* Test :
      - pose_i1_m1, pose_i1_m2
      - s_m1, alpha_m1, beta_m1
      - s_m2, alpha_m2, beta_m2

      Pose_i1_m1 = 1/s_m1 * alpha_m1^-1 * (pose_i1_m1 - beta_m1)
      Pose_i1_m2 = 1/s_m2 * alpha_m2^-1 * (pose_i1_m2 - beta_m2)

      diff = Pose_i1_m1 - Pose_i1_m2

    */
/*

---
    2 0003.png 0004.png
    alpha_m1
      -0.943879 0.0336532 0.328572
      0.0452759 0.998588 0.0277849
      -0.327173 0.041102 -0.94407
    beta_m1
      1.45337 0.782566 4.01418
    s_m1
      2.50449
---
    3 0002.png 0003.png 0004.png
    alpha_m2
      -0.835362 -0.0244854 0.549154
      -0.00363671 0.999232 0.0390211
      -0.549688 0.0305996 -0.834809
    beta_m2
      0.510063 0.616558 2.93794
    s_m2
      1.52224
---

GLOBAL POSES
    0002.png
      -0.834751 -0.00355065 -0.550616
      -0.0240994 0.999257 0.0300918
      0.5501 0.0383887 -0.834216

      1.34453 -0.452846 1.41262

    0003.png
      -0.943896 0.0443111 -0.327256
      0.0325787 0.998618 0.0412491
      0.328632 0.0282732 -0.944035

      1.05797 -0.397735 1.31378

    0004.png
      -0.978396 0.0762443 -0.192168
      0.0635704 0.995428 0.0712851
      0.196724 0.0575289 -0.97877

      0.675793 -0.394214 1.42934

---

RELATIVE POSES (pose_i_j)

    0003.png 0004.png
      0.990304 -0.0184962 -0.137682
      0.0236068 0.999088 0.0355792
      0.136898 -0.0384845 0.989837

      0.998794 -0.0280127 0.0403243


    3 0002.png 0003.png  0004.png
      0.968305 -0.0471622 -0.245278
      0.0500244 0.998733 0.00544876
      0.24471 -0.0175459 0.969438

      0.294887 0.0851058 0.36799

      0.923964 -0.055705 -0.378401
      0.073779 0.996714 0.0334227
      0.375296 -0.0587994 0.925038

      0.841244 0.100972 0.531143

*/

    // similarity trafo
    Mat3d alpha_m1;
    alpha_m1 <<  -0.943879, 0.0336532, 0.328572,
                 0.0452759, 0.998588, 0.0277849,
                 -0.327173, 0.041102, -0.94407;
    Mat3d alpha_m1_inv = alpha_m1.inverse();
    Vec3d beta_m1;
    beta_m1 << 1.45337, 0.782566, 4.01418;
    double s_m1 = 2.50449;

    Mat3d alpha_m2;
    alpha_m2 <<  -0.835362, -0.0244854, 0.549154,
                 -0.00363671, 0.999232, 0.0390211,
                -0.549688, 0.0305996, -0.834809;
    Mat3d alpha_m2_inv = alpha_m2.inverse();
    Vec3d beta_m2;
    beta_m2 << 0.510063, 0.616558, 2.93794;
    double s_m2 = 1.52224;

    //relative poses
    //motion1
    Mat3d r3_m1 = Mat3d::Identity();
    Vec3d t3_m1 = Vec3d::Zero();

    Mat3d r4_m1;
    r4_m1 <<  0.990304, -0.0184962, -0.137682,
              0.0236068, 0.999088, 0.0355792,
              0.136898, -0.0384845, 0.989837;
    Vec3d t4_m1;
    t4_m1 <<  0.998794, -0.0280127, 0.0403243;

    //motion2
    Mat3d r2_m2 = Mat3d::Identity();
    Vec3d t2_m2 = Vec3d::Zero();

    Mat3d r3_m2;
    r3_m2 <<  0.968305, -0.0471622, -0.245278,
              0.0500244, 0.998733, 0.00544876,
              0.24471, -0.0175459, 0.969438;
    Vec3d t3_m2;
    t3_m2 <<  0.294887, 0.0851058, 0.36799;

    Mat3d r4_m2;
    r4_m2 << 0.923964, -0.055705, -0.378401,
             0.073779, 0.996714, 0.0334227,
             0.375296, -0.0587994, 0.925038;
    Vec3d t4_m2;
    t4_m2 <<  0.841244, 0.100972, 0.531143;


    //absolute poses
    //Pose_i1_m1 = 1/s_m1 * alpha_m1^-1 * (pose_i1_m1 - beta_m1)
    //Pose_i1_m2 = 1/s_m2 * alpha_m2^-1 * (pose_i1_m2 - beta_m2)
    Mat3d R3_m1 = r3_m1.inverse() * alpha_m1;
    Vec3d T3_m1 = 1/s_m1 * alpha_m1_inv * (t3_m1 - beta_m1);
    cPoseBasic P3_m1(R3_m1,T3_m1,"0003.png");
    P3_m1.Show();

    Mat3d R3_m2 = r3_m2.inverse() * alpha_m2;
    Vec3d T3_m2 = 1/s_m2 * alpha_m2_inv * (t3_m2 - beta_m2);
    cPoseBasic P3_m2(R3_m2,T3_m2,"0003.png");
    P3_m2.Show();

    Mat3d R2_m2 = r2_m2.inverse() * alpha_m2;
    Vec3d T2_m2 = 1/s_m2 * alpha_m2_inv * (t2_m2 - beta_m2);
    cPoseBasic P2_m2(R2_m2,T2_m2,"0002.png");
    P2_m2.Show();


}

void gen_toy_problem_triplet(cPoseBasic* P1,cPoseBasic* P2,cPoseBasic* P3)
{

  /* input relative orientations */
  Mat3d r1;
        r1 << 1, 0, 0,
              0, 1, 0,
              0, 0, 1;
  Vec3d t1;
        t1 << 0, 0, 0;
  P1->R() = r1;
  P1->C_() = t1;

  Mat3d r2;
        r2 << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;
  Vec3d t2;
        t2 << 100, 101, 102;
  P2->R() = r2;
  P2->C_() = t2;

  Mat3d r3;
        r3 << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;
  Vec3d t3;
        t3 << 100, 101, 102;
  P3->R() = r3;
  P3->C_() = t3;

  /* input H_i, g_i */
  Mat6d H_1;
        H_1 << 1, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 1,
              0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 1;
  Vec6d g_1;
        g_1 << 1, 0, 0, 0, 0, 0;
  Mat6d H_2;
        H_2 << 1, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 1,
              0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 1;
  Vec6d g_2;
        g_2 << 1, 0, 0, 0, 0, 0;
  Mat6d H_3;
        H_3 << 1, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 1,
              0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 1;
  Vec6d g_3;
        g_3 << 1, 0, 0, 0, 0, 0;


  /* Fill-in H and g */
  Mat18d H = Mat18d::Zero();
  Vec18d g = Vec18d::Zero();

}

int main_check_similarity_trafo(int argc,char** argv)
{
    test_similarity_trafo();

    return true;
}

/* std::map<std::string,cPoseBasic*> aGlobalPoses;

//aNViewMap["0003.png-0004.png"] = motion_34;
//aNViewMap["0002.png-0003.png-0004.png"] = motion_234;

Mat3d R2;
R2 << -0.834751, -0.00355065, -0.550616,
      -0.0240994, 0.999257, 0.0300918,
       0.5501, 0.0383887, -0.834216;
Vec3d T2;
T2 << 1.34453, -0.452846, 1.41262;
cPoseBasic* P2 = new cPoseBasic(R2,T2,"0002.png");

Mat3d R3;
R3 << -0.943896, 0.0443111, -0.327256,
       0.0325787, 0.998618, 0.0412491,
       0.328632, 0.0282732, -0.944035;
Vec3d T3;
T3 <<  1.05797, -0.397735, 1.31378;
cPoseBasic* P3 = new cPoseBasic(R3,T3,"0003.png");

Mat3d R4;
R4 << -0.978396, 0.0762443, -0.192168,
       0.0635704, 0.995428, 0.0712851,
       0.196724, 0.0575289, -0.97877;
Vec3d T4;
T4 << 0.675793, -0.394214, 1.42934;
cPoseBasic* P4 = new cPoseBasic(R4,T4,"0004.png");

aGlobalPoses["0002.png"] = P2;
aGlobalPoses["0003.png"] = P3;
aGlobalPoses["0004.png"] = P4;

std::cout << aGlobalPoses.size() << "\n"; */
int main_check_propagation(int argc,char** argv)
{

  // Generate toy dataset
  cNviewPoseT<Mat6d,Vec6d>* motion_34 =0;
  cNviewPoseT<Mat12d,Vec12d>* motion_234 =0;
  test_gen_2motions(motion_34,motion_234);


  /* Propagate covariance to global frame  */
  //first motion
  test_propagate_cov(motion_34);
  //second motion
  test_propagate_cov(motion_234);


  /* Accumulate H and g */
  Mat18d global_H = Mat18d::Zero();
  Vec18d global_g = Vec18d::Zero();

  global_H.block(12,12,6,6) += motion_34->Hg_().H_();
  global_H.block(6,6,12,12) += motion_234->Hg_().H_();
  std::cout << "\n H \n" << global_H;

  global_g.block(12,0,6,1) += motion_34->Hg_().g_();
  global_g.block(6,0,12,1) += motion_234->Hg_().g_();
  std::cout << "\n g \n" << global_g  << "\n";

  /* Solve for delta
     delta =  H^-1 @ g    18 x 1
 */
 Vec18d delta = global_H.inverse() * global_g;
 std::cout << "delta=\n" << delta << "\n";


  /* Update
    Pose+ = exp(delta)^T * Pose
    [R+,T+] = exp(delta)^T * [R | T]

    T_delta = from_Rt(so3exp_map(aa), t)
    T = T_delta @ T

    so3exp_map

  */


  return true;

}



int main(int argc, char** argv)
{
    //main_check_similarity_trafo(argc,argv);

    main_check_propagation(argc,argv);

    return true;
}
