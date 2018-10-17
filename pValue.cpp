/*
  模拟椭圆的投影值生成过程
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#define PI 3.14159265354

using namespace std;

//模拟平行投影束情况下均匀椭圆的投影值
void simulation_p(double[][] &pValue,double Attenuation,double a,double b){
	double DisInterval=1;//平行束直线间隔
	double LongAxis=max(a,b);//判断椭圆长轴，决定扫描的范围
	
	for(double Theta=0;Theta<=180;Theta++){	
		for(double Distance=-LongAxis;Distance<=LongAxis;Distance+=DisInterval){
			double r2=pow(a*cos(PI*Theta/180),2)+pow(b*sin(PI*Theta/180),2);
			if(r2-Distance*Distance>=0){
				double PQ=2*a*b*sqrt(r2-Distance*Distance)/r2;
				pValue[Distance+LongAxis][Theta]=PQ*Attenuation;
			}
			else
				pValue[Distance+LongAxis][Theta]=0;
		}
	}
}

