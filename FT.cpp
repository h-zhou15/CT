#includes <iostream>
#include <cmath>
#include <vector>
#include <complex>
#define PI 3.14159265354

using namespace std;

complex<double> complex_t;

complex_t W(int N,int n,int k)
{
	return complex_t(cos(2*PI*n*k/N),-sin(2*PI*n*k/N));
}

//DFT算法
void DFT(vector<double> &x_n,vector<complex_t> &X_k)
{
	int N=x_n.size();
	X_k.clear();
	for(int k=0;k<N;k++)
	{
		complex_t t(0,0);
		for(int n=0;n<N;n++)
		{
			t+=x_n[n]*W(N,n,k);
		}
		X_k.push_back(t);
	}
}

//IDFT算法
void IDFT(vector<complex_t> &X_k,vector<double> &x_n)
{
	x_n.clear();
	int N=X_k.size();
	for(int i=0;i<n;i++)
	{
		complex_t t(0,0);
		for(int k=0;k<N;k++)
		{
			t+=X_k[k]*W(N,-i,k);
		}
		x_n.push_back(t.raal()/N);
	}
}

	
//保证N是2的n次幂
int bitlen(int N){
	int n=0;
	while((N&1)==0){
		n++;
		N>>=1;
	}
	return n;
}
 
 
int reverse_bit(int n,int len){//bit反转 
	int tmp=0;
	while(len--){
		tmp+=((n&1)<<len);
		n>>=1;
	}
	return tmp;
 
}
 
//序数重排 
void resort(vector<complex_t> &x_n,int N){
	vector<complex_t> v(x_n);
	int len=bitlen(N);
	for(int i=0;i<N;++i){
		x_n[i]=v[reverse_bit(i,len)];
	}
}
 
 
//基2,FFT算法实现,O(nlogn)的复杂度
void FFT(vector<complex_t> &x_n){
	int N=x_n.size();
	int r=bitlen(N);
	vector<complex_t> W(N);
 
	//预先计算旋转因子 
	for(int i=0;i<N;++i){
		double angle=-i*2*PI/N;
		W[i]=complex_t(cos(angle),sin(angle));
	}
	
	
	for(int k=0;k<r;++k){//迭代次数 
		for(int j=0;j<(1<<k);++j){
			int butterfly=1<<(r-k);
			int p=j*butterfly;
			int s=p+butterfly/2;
			for(int i=0;i<butterfly/2;++i){
				complex_t c=x_n[i+p]+x_n[i+s];
				x_n[i+s]=(x_n[i+p]-x_n[i+s])*W[i*(1<<k)];
				x_n[i+p]=c;
			}
		}
	}
	
	//次序重排 
	resort(x_n,N);
	for(int i=0;i<N;++i){
		cout<<format(x_n[i])<<endl;
	}
	
}
