/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Map.h"
#include "MapPoint.h"
#include "KeyFrame.h"
#include "LoopClosing.h"
#include "Frame.h"

#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"
#define Pi 3.14159265359
namespace ORB_SLAM2
{

class LoopClosing;

class Optimizer
{
public:
    void static BundleAdjustment(const std::vector<KeyFrame*> &vpKF, const std::vector<MapPoint*> &vpMP,
                                 int nIterations = 5, bool *pbStopFlag=NULL, const unsigned long nLoopKF=0,
                                 const bool bRobust = true);
    void static BundleAdjustment2(Map* pMap,
                                  int nIterations = 5, bool *pbStopFlag=NULL, const unsigned long nLoopKF=0,
                                  const bool bRobust = true);
    void static GlobalBundleAdjustemnt(Map* pMap, int nIterations=5, bool *pbStopFlag=NULL,
                                       const unsigned long nLoopKF=0, const bool bRobust = true);
    void static GlobalBundleAdjustemnt2(Map* pMap, int nIterations=5, bool *pbStopFlag=NULL,
                                       const unsigned long nLoopKF=0, const bool bRobust = true);
    void static LocalBundleAdjustment(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);
    void static LocalBundleAdjustment2(KeyFrame *pKF, bool *pbStopFlag, Map *pMap);
    int static PoseOptimization(Frame* pFrame);

    // if bFixScale is true, 6DoF optimization (stereo,rgbd), 7DoF otherwise (mono)
    void static OptimizeEssentialGraph(Map* pMap, KeyFrame* pLoopKF, KeyFrame* pCurKF,
                                       const LoopClosing::KeyFrameAndPose &NonCorrectedSim3,
                                       const LoopClosing::KeyFrameAndPose &CorrectedSim3,
                                       const map<KeyFrame *, set<KeyFrame *> > &LoopConnections,
                                       const bool &bFixScale);

    void static OptimizeEssentialGraph2(Map* pMap, KeyFrame* pLoopKF, KeyFrame* pCurKF,
                                        const LoopClosing::KeyFrameAndPose &NonCorrectedSim3,
                                        const LoopClosing::KeyFrameAndPose &CorrectedSim3,
                                        const LoopClosing::CylinderAndPose &NonCorrectedCy, const LoopClosing::CylinderAndPose &CorrectedCy,
                                        const map<KeyFrame *, set<KeyFrame *> > &LoopConnections,
                                        const bool &bFixScale);

    void static IsodiametricOptimization(Map* pMap, MapCylinder* pLastCy);
    // if bFixScale is true, optimize SE3 (stereo,rgbd), Sim3 otherwise (mono)
    static int OptimizeSim3(KeyFrame* pKF1, KeyFrame* pKF2, std::vector<MapPoint *> &vpMatches1,
                            g2o::Sim3 &g2oS12, const float th2, const bool bFixScale);
};

class CylinderFittingVertex : public g2o::BaseVertex<5, Eigen::Matrix<double,5,1>> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual void setToOriginImpl() override {
        _estimate << 1.0, 2.0, 2.0, 3.0, 9.0;
    }


    virtual void oplusImpl(const double *update) override {
        Eigen::Matrix<double,5,1> last=_estimate;
        _estimate += Eigen::Matrix<double,5,1>(update);
        if(floor(_estimate[0]*2/Pi)>floor(last[0]*2/Pi)) _estimate[0]=floor(last[0]*2/Pi+1)*Pi/2;
        if(floor(_estimate[1]*2/Pi)>floor(last[1]*2/Pi)) _estimate[1]=floor(last[1]*2/Pi+1)*Pi/2;
    }

    virtual bool read(std::istream &in) {}

    virtual bool write(std::ostream &out) const {}
};


class CylinderFittingEdge : public g2o::BaseBinaryEdge<1, double, g2o::VertexSBAPointXYZ, CylinderFittingVertex> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual void computeError() override {
        const g2o::VertexSBAPointXYZ *vPoint = static_cast<const g2o::VertexSBAPointXYZ *> (_vertices[0]);
        const CylinderFittingVertex *v = static_cast<const CylinderFittingVertex *> (_vertices[1]);
        const Eigen::Matrix<double,5,1> abc = v->estimate();
        double a = abc[0];
        double b = abc[1];
        double Px = abc[2];
        double Py = abc[3];
        Eigen::Vector3d Lcy(cos(a)*sin(b),sin(a)*sin(b), cos(b));
        Eigen::Vector3d Ocy(cos(a)*cos(b)*Px-sin(a)*Py,sin(a)*cos(b)*Px+cos(a)*Py,-sin(b)*Px);
        _error(0, 0) =(Lcy.cross(vPoint->estimate()+Ocy)).norm()-abc[4];
    }

    virtual void linearizeOplus() override {
        const g2o::VertexSBAPointXYZ *vPoint = static_cast<const g2o::VertexSBAPointXYZ *> (_vertices[0]);
        const CylinderFittingVertex *v = static_cast<const CylinderFittingVertex *> (_vertices[1]);
        const Eigen::Matrix<double,5,1> abc = v->estimate();
        double a = abc[0];
        double b = abc[1];
        double Px = abc[2];
        double Py = abc[3];
        double X=vPoint->estimate()[0];
        double X2=pow(X,2);
        double Y=vPoint->estimate()[1];
        double Y2=pow(Y,2);
        double Z=vPoint->estimate()[2];
        double Z2=pow(Z,2);
        double ca=cos(a);
        double ca2=pow(cos(a),2);
        double sa=sin(a);
        double sa2=pow(sin(a),2);
        double cb=cos(b);
        double cb2=pow(cos(b),2);
        double sb=sin(b);
        double sb2=pow(sin(b),2);

        double delta=1/sqrt(cb2*(X2+Y2)+pow(Px,2)+pow(Py,2)+sb2*Z2-2*sin(b)*Px*Z
                            +2*sa*cb*Px*Y+2*ca*cb2*Py*Y-2*sa*cb*sb*Y*Z
                            +2*ca*cb*Px*X-2*sa*cb2*Py*X-2*ca*cb*sb*X*Z
                            +sa2*sb2*X2+ca2*sb2*Y2-2*sa*sb2*Py*X-2*ca*sa*sb2*X*Y
                            +2*cos(a)*pow(sin(b),2)*Py*Y);
        _jacobianOplusXi[0] = cb2*X+ca*cb*Px-sa*cb2*Py-ca*cb*sb*Z
                              +sa2*sb2*X-sa*sb2*Py-ca*sa*sb2*Y;
        _jacobianOplusXi[0]=_jacobianOplusXi[0]*delta;
        _jacobianOplusXi[1] = cb2*Y+sa*cb*Px+ca*Py-sa*cb*sb*Z
                              +ca2*sb2*Y-ca*sa*sb2*X;
        _jacobianOplusXi[1]=_jacobianOplusXi[1]*delta;
        _jacobianOplusXi[2] = sb2*Z-sb*Px-sa*cb*sb*Y-ca*cb*sb*X;
        _jacobianOplusXi[2]=_jacobianOplusXi[2]*delta;


        _jacobianOplusXj[0] =ca*cb*Px*Y-sa*cb2*Py*Y-ca*cb*sb*Y*Z
                             -sa*cb*Px*X-ca*cb2*Py*X+sa*cb*sb*X*Z
                             +sa*ca*sb2*X2-ca*sa*sb2*Y2-ca*sb2*Py*X-cos(2*a)*sb2*X*Y-sa*sb2*Py*Y;
        _jacobianOplusXj[0]=_jacobianOplusXj[0]*delta;
        _jacobianOplusXj[1] = -cb*sb*(X2+Y2)+cb*sb*Z2-cb*Px*Z
                              -sa*sb*Px*Y-2*ca*cb*sb*Py*Y-sa*cos(2*b)*Y*Z
                              -ca*sb*Px*X+2*sa*cb*sb*Py*X-ca*cos(2*b)*X*Z
                              +sa2*cb*sb*X2+ca2*cb*sb*Y2-2*sa*cb*sb*Py*X
                              -2*ca*sb*cb*sb*X*Y+2*ca*cb*sb*Py*Y;
        _jacobianOplusXj[1]=_jacobianOplusXj[1]*delta;
        _jacobianOplusXj[2] = Px-sb*Z+sa*cb*Y+ca*cb*X;
        _jacobianOplusXj[2]=_jacobianOplusXj[2]*delta;
        _jacobianOplusXj[3] = Py+ca*Y-sa*X;
        _jacobianOplusXj[3]=_jacobianOplusXj[3]*delta;
        _jacobianOplusXj[4] = -1;
    }

    virtual bool read(std::istream &in) {}

    virtual bool write(std::ostream &out) const {}
};

} //namespace ORB_SLAM

#endif // OPTIMIZER_H
