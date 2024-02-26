//
// Created by uos on 01/12/2021.
//

#ifndef ORB_SLAM2_MAPCYLINDER_H
#define ORB_SLAM2_MAPCYLINDER_H

#include"Map.h"

#include<opencv2/core/core.hpp>

namespace ORB_SLAM2
{
class Map;
//class MapPoint;
class MapCylinder{
public:
    MapCylinder(const cv::Mat &Para, Map* pMap);
    MapCylinder(const cv::Mat &Para);
    MapCylinder(Map *pMap);


    cv::Mat GetWPara();
    void SetWorldPara(const cv::Mat &Para);
    void SetWorldMat(const cv::Mat &TwcyH);
    void SetCyMat(const cv::Mat &TcywH);
    void AddMapPoint(MapPoint* pMP);
    void UpdateMapCylinder(const cv::Mat &Para);
    void CalculateLength();
    void CalculateLength(list<MapPoint*> lCyPoint, cv::Mat curKFpos);
    void ToCyMatrix();
    void GetCyWMat(cv::Mat &Rwcy, cv::Mat &Pwcy, cv::Mat &Twcy);
    void GetCyCyMat(cv::Mat &Rcyw, cv::Mat &Pcyw, cv::Mat &Tcyw);
    cv::Mat GetPose();
    cv::Mat GetRotation();
    cv::Mat GetTranslation();
    cv::Mat GetPoseInverse();
    cv::Mat GetDirection();
    float GetRadius();
    float GetStart();
    float GetEnd();
    void updateLength(float position);
    void updateLength(float position1, float position2);
    list<KeyFrame*> GetCylindricalKF();
    void AddCylindricalKF(KeyFrame* pKF);
    void UpdateCylindricalKFs();
    void CorrectScale(double s);
    std::map<KeyFrame*,size_t> GetObservations();
    void cyPreparation(std::list<MapPoint *> lFixedMapPoints);
    void ransacGroup(std::vector<std::vector<size_t> > &mvSets, int mMaxIterations, size_t N);
    void initialValue( list<KeyFrame*> lLocalKeyFrames);
    void SetBadFlag();
public:
    long unsigned int mnId;
    static long unsigned int nNextId;
    long unsigned int mnBALocalForKF;
    long unsigned int mnLoopClosingForKF;
    double inverseVariance;
    list<MapPoint*> mlpMapPoints;
    bool bActive;               // start to estimate or help triangulate
    bool bBuilt;
    list<float> mlMPsOnAxis;     //Projections of points on the axis
    cv::Mat mParaGBA;
    long unsigned int mnBAGlobalForCy;

    float fx;
    float fy;
    float cx;
    float cy;
    float invfx;
    float invfy;

    std::vector<cv::Mat> vFixedMapPointPos;
    std::vector<cv::Mat> vKFpos1;
    std::vector<cv::Mat> vKFpos2;
    std::vector<cv::KeyPoint> vMpobs1;
    std::vector<cv::KeyPoint> vMpobs2;
    std::vector<bool> inlierList;
    std::vector<float> vSigmaSquare2;
protected:
    cv::Mat mPara;     //5 parameters, world coordinate
    cv::Mat mTwcy;     //T metrix, world coordinate
    cv::Mat mRwcy;      //Rotation matrix, world coordinate
    cv::Mat mPwcy;       //position, world coordinate
    cv::Mat mTcyw;     //T metrix, world coordinate
    cv::Mat mRcyw;      //Rotation matrix, world coordinate
    cv::Mat mPcyw;       //position, world coordinate
    cv::Mat mDir;           //direction, world coordinate

    float mfstart;
    float mfend;
    float mMaxMPonAxis;
    float mMinMPonAxis;
    list<KeyFrame*> mlpCylindricalKF;
    std::map<KeyFrame*,size_t> mObservations;
    Map* mpMap;
};

} //namespace ORB_SLAM

#endif //ORB_SLAM2_MAPCYLINDER_H
