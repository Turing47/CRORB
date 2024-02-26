//
// Created by uos on 01/12/2021.
//

#include "MapCylinder.h"

namespace ORB_SLAM2 {

long unsigned int MapCylinder::nNextId=0;

MapCylinder::MapCylinder(const cv::Mat &Para, Map* pMap):mnBALocalForKF(0),mnLoopClosingForKF(0),mpMap(pMap),bActive(false),bBuilt(false),mfstart(0),mfend(0),mMaxMPonAxis(-10000000),mMinMPonAxis(10000000){
    Para.copyTo(mPara);
    ToCyMatrix();
    mnId=nNextId++;
    inverseVariance=pow(39.2/mPara.at<float>(4),2);
}

MapCylinder::MapCylinder(const cv::Mat &Para):bActive(false),bBuilt(false),mfstart(0),mfend(0),mMaxMPonAxis(-10000000),mMinMPonAxis(10000000),mnBALocalForKF(0),mnLoopClosingForKF(0){
    mPara=Para;
    ToCyMatrix();
    mnId=nNextId++;
        inverseVariance=pow(39.2/mPara.at<float>(4),2);
}

MapCylinder::MapCylinder(Map* pMap):mnBALocalForKF(0),mnLoopClosingForKF(0),mpMap(pMap),bActive(false),bBuilt(false),mfstart(0),mfend(0),mMaxMPonAxis(-10000000),mMinMPonAxis(10000000){
    mPara=cv::Mat::ones(5,1,CV_32F);
    mnId=nNextId++;
}

cv::Mat MapCylinder::GetWPara(){
    return mPara.clone();
}
void MapCylinder::SetWorldPara(const cv::Mat &Para){
    Para.copyTo(mPara);
    ToCyMatrix();
    inverseVariance=pow(39.2/mPara.at<float>(4),2);
}

void MapCylinder::SetWorldMat(const cv::Mat &TwcyH) {
    cv::Mat wDir=TwcyH.rowRange(0,3).col(2);
    wDir.copyTo(mDir);
    float l=wDir.at<float>(0);
    float m=wDir.at<float>(1);
    float n=wDir.at<float>(2);
    cv::Mat PwcyH=TwcyH.rowRange(0,3).col(3);
    float x=PwcyH.at<float>(0);
    float y=PwcyH.at<float>(1);
    float z=PwcyH.at<float>(2);
    float t=-l*x-m*y-n*z;
    cv::Mat Pwcy(3,1,CV_32F);

    Pwcy.at<float>(0)=l*t+x;
    Pwcy.at<float>(1)=m*t+y;
    Pwcy.at<float>(2)=n*t+z;

    float devOnAxis=-abs(t)/t*sqrt(pow((Pwcy.at<float>(0)-x),2)+pow((Pwcy.at<float>(1)-y),2)+pow((Pwcy.at<float>(2)-z),2));
    mfstart+=devOnAxis;
    mfend+=devOnAxis;
    for(std::list<float>::iterator lit=mlMPsOnAxis.begin(), lend=mlMPsOnAxis.end(); lit!=lend; lit++){
        *lit+=devOnAxis;
    }

    float a,b,ratio;
    b=acos(n);
    ratio=l/sin(b);
    if(ratio>0.999999) a=0;
    else if(ratio<-0.999999) a=M_PI;
    else a=acos(ratio);
    if(m<0) a=-a;

    float ca=cos(a);
    float sa=sin(a);
    float cb=cos(b);
    float sb=sin(b);
    cv::Mat Rwcy(3,3,CV_32F);
    Rwcy.at<float>(0,0)=cb*ca;
    Rwcy.at<float>(0,1)=-sa;
    Rwcy.at<float>(0,2)=sb*ca;

    Rwcy.at<float>(1,0)=cb*sa;
    Rwcy.at<float>(1,1)=ca;
    Rwcy.at<float>(1,2)=sb*sa;

    Rwcy.at<float>(2,0)=-sb;
    Rwcy.at<float>(2,1)=0;
    Rwcy.at<float>(2,2)=cb;
    Rwcy.copyTo(mRwcy);
    cv::Mat Twcy=cv::Mat::eye(4,4,CV_32F);
    Rwcy.copyTo(Twcy.rowRange(0,3).colRange(0,3));
    Pwcy.copyTo(Twcy.rowRange(0,3).col(3));
    Pwcy.copyTo(mPwcy);
    Twcy.copyTo(mTwcy);

    cv::Mat Tcyw=cv::Mat::eye(4,4,CV_32F);
    cv::Mat Rcyw=Rwcy.t();
    Rcyw.copyTo(mRcyw);
    cv::Mat Pcyw=-Rcyw*Pwcy;
    Rcyw.copyTo(Tcyw.rowRange(0,3).colRange(0,3));
    Pcyw.copyTo(Tcyw.rowRange(0,3).col(3));
    Tcyw.copyTo(mTcyw);
    mPara.at<float>(0)=a;
    mPara.at<float>(1)=b;
    mPara.at<float>(2)=Pcyw.at<float>(0);
    mPara.at<float>(2)=Pcyw.at<float>(1);
}

void MapCylinder::SetCyMat(const cv::Mat &TcywH){
    cv::Mat wDir=TcywH.row(2).colRange(0,3);
    float l=wDir.at<float>(0);
    float m=wDir.at<float>(1);
    float n=wDir.at<float>(2);
    cv::Mat PwcyH=-TcywH.rowRange(0,3).colRange(0,3).t()*TcywH.rowRange(0,3).col(3);

    float a,b,ratio;
    b=acos(n);
    ratio=l/sin(b);
    if(ratio>0.999999) a=0;
    else if(ratio<-0.999999) a=M_PI;
    else a=acos(ratio);
    if(m<0) a=-a;

    float ca=cos(a);
    float sa=sin(a);
    float cb=cos(b);
    float sb=sin(b);
    cv::Mat Rwcy(3,3,CV_32F);
    Rwcy.at<float>(0,0)=cb*ca;
    Rwcy.at<float>(0,1)=-sa;
    Rwcy.at<float>(0,2)=sb*ca;

    Rwcy.at<float>(1,0)=cb*sa;
    Rwcy.at<float>(1,1)=ca;
    Rwcy.at<float>(1,2)=sb*sa;

    Rwcy.at<float>(2,0)=-sb;
    Rwcy.at<float>(2,1)=0;
    Rwcy.at<float>(2,2)=cb;
    Rwcy.copyTo(mRwcy);

    cv::Mat Pcyw(3,1,CV_32F);
    Pcyw=-Rwcy.t()*PwcyH;
    float devOnAxis=-Pcyw.at<float>(2,0);
    Pcyw.at<float>(2,0)=0;
    Pcyw.copyTo(mPcyw);

    mfstart+=devOnAxis;
    mfend+=devOnAxis;
    for(std::list<float>::iterator lit=mlMPsOnAxis.begin(), lend=mlMPsOnAxis.end(); lit!=lend; lit++){
        *lit+=devOnAxis;
    }

    cv::Mat Pwcy(3,1,CV_32F);
    Pwcy=-Rwcy*Pcyw;
    cv::Mat Twcy=cv::Mat::eye(4,4,CV_32F);
    Rwcy.copyTo(Twcy.rowRange(0,3).colRange(0,3));
    Pwcy.copyTo(Twcy.rowRange(0,3).col(3));
    Pwcy.copyTo(mPwcy);
    Twcy.copyTo(mTwcy);

    cv::Mat Tcyw=cv::Mat::eye(4,4,CV_32F);
    cv::Mat Rcyw=Rwcy.t();
    Rcyw.copyTo(mRcyw);
    Rcyw.copyTo(Tcyw.rowRange(0,3).colRange(0,3));
    Pcyw.copyTo(Tcyw.rowRange(0,3).col(3));
    Tcyw.copyTo(mTcyw);
    mPara.at<float>(0)=a;
    mPara.at<float>(1)=b;
    mPara.at<float>(2)=Pcyw.at<float>(0);
    mPara.at<float>(3)=Pcyw.at<float>(1);
}

void MapCylinder::CorrectScale(double s) {
    mPara.at<float>(4)*=(1./s);
    mfstart*=(1./s);
    mfend*=(1./s);
    for(std::list<float>::iterator lit=mlMPsOnAxis.begin(), lend=mlMPsOnAxis.end(); lit!=lend; lit++){
        *lit*=(1./s);
    }
}

void MapCylinder::AddMapPoint(MapPoint* pMP){
    mlpMapPoints.push_back(pMP);
    if(bBuilt==1){
        float dLength=pMP->GetWorldPos().dot(mDir);
        mlMPsOnAxis.push_back(dLength);
        if(dLength>mMaxMPonAxis) mMaxMPonAxis=dLength;
        else if(dLength<mMinMPonAxis) mMinMPonAxis=dLength;

        KeyFrame* pRefKF=pMP->GetReferenceKeyFrame();
        mObservations[pRefKF]++;
        if(mObservations[pRefKF]==200) AddCylindricalKF(pRefKF);
    }
}

void MapCylinder::UpdateMapCylinder(const cv::Mat &Para){
    Para.copyTo(mPara);
    ToCyMatrix();
    inverseVariance=pow(39.2/mPara.at<float>(4),2);
//    inverseVariance=pow(784/39/mPara.at<float>(4),2);
}

void MapCylinder::CalculateLength(){
    size_t n=39;
    vector<int> histogramMPonAxis(n+1,0);
    float r=abs(mPara.at<float>(4));

    for(std::list<float>::iterator lit=mlMPsOnAxis.begin(), lend=mlMPsOnAxis.end(); lit!=lend; lit++) {
        float MPonAxis = *lit;
        size_t index = floor((MPonAxis-mMinMPonAxis) / r);
        if(index>n){
            n=index+1;
            histogramMPonAxis.resize(n+40);
        }
        histogramMPonAxis[index]++;
    }

    int thresholdAcc=ceil(mlMPsOnAxis.size()/10);
    for(size_t i=0;i<histogramMPonAxis.size();i++){
        if (histogramMPonAxis[i]>thresholdAcc) {
            mfstart=(i+0.5)*r+mMinMPonAxis;
            break;
        }
    }

    auto maxIndex= max_element(histogramMPonAxis.begin(),histogramMPonAxis.end());
    int sum=accumulate(histogramMPonAxis.begin(),maxIndex,0);
    unsigned int thresholdSingle=thresholdAcc*2;
    thresholdAcc*=5;

    float averageNum;
    for(int i=maxIndex-histogramMPonAxis.begin();i<histogramMPonAxis.size();i++) {
        sum += histogramMPonAxis[i];
        averageNum=0;
        if (sum < thresholdAcc) {
            for (int j = -1; j < 2; j++) {
                if (i + j < 0) averageNum += histogramMPonAxis[0];
                else if (i + j >= histogramMPonAxis.size())
                    averageNum += histogramMPonAxis[histogramMPonAxis.size() - 1];
                else averageNum += histogramMPonAxis[i + j];
            }
            if (averageNum > thresholdSingle) {
                mfend = mMinMPonAxis + (i + 0.5) * r;
                break;
            }
        }
        else{
            mfend = mMinMPonAxis + (i + 0.5) * r;
            break;
        }
    }
}

void MapCylinder::CalculateLength(list<MapPoint*> lCyPoint, cv::Mat curKFpos){
    float fMPonAxis;
    float fminMPonAxis=curKFpos.dot(mDir);
    list<float> LocalMPsOnAxis;
    float distant;
    bool bdistant=0;
    vector<vector<float>> histogramMPonAxis(40);
    for(list<MapPoint*>::iterator lit=lCyPoint.begin(), lend=lCyPoint.end(); lit!=lend; lit++) {
        MapPoint* pMP = *lit;
        pMP->SetCylinder(this);
        mlpMapPoints.push_back(pMP);

        fMPonAxis=pMP->GetWorldPos().dot(mDir);
        mlMPsOnAxis.push_back(fMPonAxis);
        LocalMPsOnAxis.push_back(fMPonAxis);

        KeyFrame* pRefKF=pMP->GetReferenceKeyFrame();
        if(mObservations.count(pRefKF)){
            mObservations[pRefKF]++;
            if(mObservations[pRefKF]==200)  AddCylindricalKF(pRefKF);
        }
        else {
            mObservations[pRefKF] = 1;
        }
    }

    float interval=abs(mPara.at<float>(4))/4;


    for(std::list<float>::iterator lit=LocalMPsOnAxis.begin(), lend=LocalMPsOnAxis.end(); lit!=lend; lit++) {
        float MPonAxis = *lit;
        size_t index = floor(abs(MPonAxis-fminMPonAxis) / interval);
        if(index>=histogramMPonAxis.size()) continue;
        histogramMPonAxis[index].push_back(MPonAxis);
    }
    std::vector<size_t>lhistogramNum(histogramMPonAxis.size());
    size_t maxNumMPinHistogram=histogramMPonAxis[0].size();
    size_t maxIndexMPinHistogram=0;
    for(size_t i=0;i<lhistogramNum.size();i++){
        lhistogramNum[i]=histogramMPonAxis[i].size();
       if(lhistogramNum[i]>maxNumMPinHistogram) {
            maxNumMPinHistogram=lhistogramNum[i];
            maxIndexMPinHistogram=i;
        }
    }
//
    int thresholdIndex=ceil(LocalMPsOnAxis.size()*0.1)>maxNumMPinHistogram?ceil(LocalMPsOnAxis.size()*0.1):maxNumMPinHistogram;
    float averageNum=0;

    int sum=accumulate(lhistogramNum.begin(),lhistogramNum.begin()+maxIndexMPinHistogram,0);

    unsigned int thresholdsum=LocalMPsOnAxis.size()*0.8;
    for(int i=maxIndexMPinHistogram;i<histogramMPonAxis.size();i++){
        sum+=int(lhistogramNum[i]);
        if (sum<thresholdsum) {
            if(i-1<0) averageNum=(lhistogramNum[0]+lhistogramNum[1])/2;
            else if(i+1>=lhistogramNum.size()) averageNum=(lhistogramNum[lhistogramNum.size()-1]+lhistogramNum[lhistogramNum.size()-2])/2;
            else averageNum+=(lhistogramNum[i-1]+lhistogramNum[i]+lhistogramNum[i+1])/3;

            if(averageNum<thresholdIndex && histogramMPonAxis[i].size()!=0) {
                distant = *max_element(histogramMPonAxis[i].begin(), histogramMPonAxis[i].end());
                if(distant-fminMPonAxis<0) distant = *min_element(histogramMPonAxis[i].begin(), histogramMPonAxis[i].end());
                bdistant=1;
                break;
            }
        }
        else {
            distant = *max_element(histogramMPonAxis[i].begin(), histogramMPonAxis[i].end());
            if(distant-fminMPonAxis<0) distant = *min_element(histogramMPonAxis[i].begin(), histogramMPonAxis[i].end());
            bdistant=1;
            break;
        }
    }

    if (bdistant==1){
        if(distant-fminMPonAxis>0 && distant > mfend) mfend = distant;
        if(distant-fminMPonAxis<0 && distant < mfstart) mfstart = distant;
    }
    return;
}

void MapCylinder::ToCyMatrix(){
    double ca=cos(mPara.at<float>(0));
    double sa=sin(mPara.at<float>(0));
    double cb=cos(mPara.at<float>(1));
    double sb=sin(mPara.at<float>(1));
    cv::Mat Rwcy(3,3,CV_32F);
    Rwcy.at<float>(0,0)=cb*ca;
    Rwcy.at<float>(0,1)=-sa;
    Rwcy.at<float>(0,2)=sb*ca;

    Rwcy.at<float>(1,0)=cb*sa;
    Rwcy.at<float>(1,1)=ca;
    Rwcy.at<float>(1,2)=sb*sa;

    Rwcy.at<float>(2,0)=-sb;
    Rwcy.at<float>(2,1)=0;
    Rwcy.at<float>(2,2)=cb;
    Rwcy.copyTo(mRwcy);
    mRcyw=Rwcy.t();

    cv::Mat Pcyw(3,1,CV_32F);
    Pcyw.at<float>(0,0)=mPara.at<float>(2);
    Pcyw.at<float>(1,0)=mPara.at<float>(3);
    Pcyw.at<float>(2,0)=0;
    mPcyw=Pcyw;
    mPwcy=-mRwcy*Pcyw;

    cv::Mat Twcy=cv::Mat::eye(4,4,CV_32F);
    mRwcy.copyTo(Twcy.rowRange(0,3).colRange(0,3));
    mPwcy.copyTo(Twcy.rowRange(0,3).col(3));
    mTwcy=Twcy;

    cv::Mat Tcyw=cv::Mat::eye(4,4,CV_32F);
    mRcyw.copyTo(Tcyw.rowRange(0,3).colRange(0,3));
    mPcyw.copyTo(Tcyw.rowRange(0,3).col(3));
    mTcyw=Tcyw;

    cv::Mat Dir(3,1,CV_32F);
    Dir.at<float>(0,0)=ca*sb;
    Dir.at<float>(1,0)=sa*sb;
    Dir.at<float>(2,0)=cb;
    mDir=Dir;
}

void MapCylinder::GetCyWMat(cv::Mat &Rwcy, cv::Mat &Pwcy, cv::Mat &Twcy) {
    Rwcy=mRwcy.clone();
    Pwcy=mPwcy.clone();
    Twcy=mTwcy.clone();
}

void MapCylinder::GetCyCyMat(cv::Mat &Rcyw, cv::Mat &Pcyw, cv::Mat &Tcyw) {
    Rcyw=mRcyw.clone();
    Pcyw=mPcyw.clone();
    Tcyw=mTcyw.clone();
}

cv::Mat MapCylinder::GetPose(){
    return mTcyw.clone();
}

cv::Mat MapCylinder::GetRotation(){
    return mRcyw.clone();
}

cv::Mat MapCylinder::GetTranslation(){
    return mPcyw.clone();
}

cv::Mat MapCylinder::GetPoseInverse() {
    return mTwcy.clone();
}

cv::Mat MapCylinder::GetDirection(){
    return mDir.clone();
}

float MapCylinder::GetRadius(){
    return mPara.at<float>(4);
}

float MapCylinder::GetStart(){
    return mfstart;
}
float MapCylinder::GetEnd(){
    return mfend;
}

void MapCylinder::updateLength(float position){
    if(position<mfstart) mfstart=position;
    else if(position>mfend) mfend=position;
}

void MapCylinder::updateLength(float position1, float position2){
    float large;
    float small;
    if(position1>=position2) {
        large=position1;
        small=position2;
    }
    else{
        large=position2;
        small=position1;
    }
    if(large>mfend) mfend=large;
    if(small<mfstart) mfstart=small;
}

list<KeyFrame*> MapCylinder::GetCylindricalKF() {
    return mlpCylindricalKF;
}

void MapCylinder::UpdateCylindricalKFs(){
    for(list<MapPoint*>::iterator lit=mlpMapPoints.begin(), lend=mlpMapPoints.end(); lit!=lend; lit++) {
        KeyFrame* pRefKF=(*lit)->GetReferenceKeyFrame();
        mObservations[pRefKF]++;
        if(mObservations[pRefKF]==200) AddCylindricalKF(pRefKF);
    }
}

std::map<KeyFrame*,size_t> MapCylinder::GetObservations(){
    return mObservations;
}

void MapCylinder::AddCylindricalKF(KeyFrame* pKF) {
    mlpCylindricalKF.push_back(pKF);
    pKF->AddMapCylinder(this);
}

void MapCylinder::initialValue( list<KeyFrame*> lLocalKeyFrames){
    int num=lLocalKeyFrames.size();
    cv::Mat data(3,num,CV_32F);
    cv::Mat twc;
    int i=0;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++){
        twc=(*lit)->GetCameraCenter();
        data.at<float>(0,i)=twc.at<float>(0,0);
        data.at<float>(1,i)=twc.at<float>(1,0);
        data.at<float>(2,i)=twc.at<float>(2,0);
        i++;
    }

    cv::PCA pca(data,cv::Mat(),1,1);

    float b=acos(pca.eigenvectors.at<float>(2));
    float ratio=pca.eigenvectors.at<float>(0)/sin(b);
    float a;
    if(ratio>0.999999) a=0;
    else if(ratio<-0.999999) a=M_PI;
     else a=acos(ratio);
    if(pca.eigenvectors.at<float>(1)<0) a=-a;

    float ca=cos(a);
    float sa=sin(a);
    float cb=cos(b);
    float sb=sin(b);
    cv::Mat Rwcy(3,3,CV_32F);
    Rwcy.at<float>(0,0)=cb*ca;
    Rwcy.at<float>(0,1)=-sa;
    Rwcy.at<float>(0,2)=sb*ca;

    Rwcy.at<float>(1,0)=cb*sa;
    Rwcy.at<float>(1,1)=ca;
    Rwcy.at<float>(1,2)=sb*sa;

    Rwcy.at<float>(2,0)=-sb;
    Rwcy.at<float>(2,1)=0;
    Rwcy.at<float>(2,2)=cb;

    cv::Mat Pcyw(3,1,CV_32F);
    Pcyw=-Rwcy.t()*pca.mean;
    mPara.at<float>(0)=a;
    mPara.at<float>(1)=b;
    mPara.at<float>(2)=Pcyw.at<float>(0);
    mPara.at<float>(3)=Pcyw.at<float>(1);
    mPara.at<float>(4)=0.3;
}

void MapCylinder::SetBadFlag(){
    for(list<MapPoint*>::iterator lit=mlpMapPoints.begin(), lend=mlpMapPoints.end(); lit!=lend; lit++){
        (*lit)->mpMapCylinder=static_cast<MapCylinder*>(NULL);
    }
    for(list<KeyFrame*>::iterator lit=mlpCylindricalKF.begin() , lend=mlpCylindricalKF.end(); lit!=lend; lit++){
        (*lit)->AddMapCylinder(NULL);
    }
    mpMap->EraseMapCylinder(this);
}
}
