using namespace std;
#include<iostream>
#include<fstream>
void normalize(double [], int);
void convolution (double arr[],double coef[],int length,int len2);
void square (double [],int );
double thresh(double [],int);
void interval(double [],double,int [],int [],int*,int*,int);
void maxima(double [],int,int,int*);
void minima(double [],int,int,int*);
double avg(double,double);
double pulserate(int,int,int);
double mod(double);
void tachy_brady(int);
void bloodpressure(double,double*,double*);
int main(){

	//////////////////////////////////////////
	/////////loading file////////////////////
	////////////////////////////////////////
	ifstream out;
	//ofstream in,inl,inr,rlc,qlc,slc;
	//ofstream plc;	
	out.open("VBIG.txt");
	double ecg[25000]={0};
	double data[25000]={0};
	int rloc[6000]={0};
	int i=0;
	int start=0,end=0;
	int vfcount=0;						//consecutive ventricular fibrillation/flutter count
	int pvccount=0;						//consecutive PVCs count
	int blockcount=0;					//consecutive heartblock beat count
	int pvcflag=0; //arrhythmia flag
	/* zero for normal beat
one for pvc beat
		two for any other beat*/

	int bicount=0;
	int bistart=0;
	int biend=0;
	int biendflag=0;   ///////identifies end of bigeminy episode
	int pos=1;
	int length=0,len2=0,len3=0;
	int llength=0,rlength=0;		//array lengths of left and right QRS complex limits
	double systole=0,dystole=0;		//blood pressure variables
	double PR=0;		//pulse rate
	int fs=360;			//sampling frequency 
	double sum=0;		//sum for calculating mean
	double mean=0;		//mean of ecg
	double thr=0;		//detection threshold
	double lp_coef[13]={0,0,1,2,3,4,5,6,5,4,3,2,1};		//low pass filter coefficients
	double hp_coef[33]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,31,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};		//high pass filer coefficients
	double h_der[5]={0.125,0.25,0,-0.25,-0.125}; 	//derivative coeficients 
	double h_int[31]={ 0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323,0.0323};		//integration coefficients
	int left[6000]={0};  //QRS complex start positions 
	int right[6000]={0};	//QRS complex end positions
	double rint[6000]={0};   //RR interval
	
	//////////reading from file////////////////////
	
	while(!out.eof()){
		out>>ecg[i];
		data[i]=ecg[i];
		i++;

	}
	length=i-1;

	///////////////////////////////////////
	//////////calculating mean////////////
	/////////////////////////////////////
	for(i=0;i<length;i++){
		sum+=ecg[i];
	}
	mean=sum/(length-1);

	////cancelling dc component///////
	for(i=0;i<length;i++){
		ecg[i]-=mean;
		data[i]-=mean;
	}
	normalize(ecg,length);	/////////////// normalize to one//////////////
	normalize(data,length);

	///////////////////////////////////////////////
	////////////////low pass filter////////////////
	//////////////////////////////////////////////
	len2=sizeof(lp_coef)/sizeof(double);
	convolution(ecg,lp_coef,length,len2);
	normalize(ecg,length+len2-1);

	/////////////////////////////////////////////
	//////////////High pass filter//////////////
	///////////////////////////////////////////
	len3=len2+length-1;
	len2=sizeof(hp_coef)/sizeof(double);	
	convolution(ecg,hp_coef,len3,len2);
	normalize(ecg,len3+len2-1);

	//////////////////////////////////////
	//////////derivative/////////////////
	////////////////////////////////////
	len3+=len2-1;
	len2=sizeof(h_der)/sizeof(double);
	convolution(ecg,h_der,len3,len2);
	normalize(ecg,length);

	///////////////////////////////////
	//////////squaring////////////////
	/////////////////////////////////
	square(ecg,length);
	normalize(ecg,length);


	////////////////////////////////////
	////////////integration////////////
	//////////////////////////////////
	len2=sizeof(h_int)/sizeof(double);
	convolution(ecg,h_int,length,len2);
	normalize(ecg,length);
	

	///////////////////////////////////////////
	///////////calculating threshold//////////
	/////////////////////////////////////////
	thr=thresh(ecg,length);

	//////////////////////////////////////////
	////////finding QRS-complex//////////////
	////////////////////////////////////////

	interval(ecg,thr,left,right,&llength,&rlength,length);
	
	///////////////////////////////////////
	/////////Q R S locations//////////////
	/////////////////////////////////////
	for(i=0;i<llength;i++){
		maxima(data,left[i],right[i],&rloc[i]);
		rloc[i]++;			///////////////offset of one from matlab results due to unknown reason in all values;
	}


	
	//////////////////////////////////////
	/////Pulse rate parameters///////////
	////////////////////////////////////
	
	PR=pulserate(length,llength,fs);
	cout<<"Pulse rate= "<<(int)PR<<endl;
	tachy_brady(PR);
	bloodpressure(PR,&systole,&dystole);
	cout<<"Systole= "<<(int)systole<<endl;
	cout<<"Dystole= "<<(int)dystole<<endl;
	
	////////////////finding RR intervals///////////
	
	for(i=0;i<llength-1;i++){
		rint[i]=((double)rloc[i+1]-(double)rloc[i])/(double)fs;
	}
	
	//////////////////////////////////////////////////////
	///////arrhythmia detection based on RR interval/////
	////////////////////////////////////////////////////

	for(i=0;i<llength-3;i++){
		vfcount=0;
		pvccount=0;
		blockcount=0;
		start=0;
		end=0;

		////////////////////////////////////////////
		/////////////////V Fib/////////////////////
		//////////////////////////////////////////

		if(rint[i+1]<0.6&&(1.8*rint[i+1])<rint[i]){
			vfcount=1;
			pvccount=0;
			pvcflag=2;
			blockcount=0;
			start=i;
			i++;
			while((rint[i]<0.7&&rint[i+1]<0.7&&rint[i+2]<0.7)||(rint[i]+rint[i+1]+rint[i+2])<1.7){
				vfcount++;
				i++;
				if(i>=llength-3){
					i=llength-3;
					end=i;
					break;
				}
				end=i;
			}
			if(vfcount>=4){
				cout<<"Ventricular Fibrillation detected from "<<(double)rloc[start]/double(fs)<<" to "<<(double)rloc[end]/(double)fs<<endl;
			}
		}

		////////////////////////////////////////////
		/////////////////PVC///////////////////////
		//////////////////////////////////////////

		else if(((1.15*rint[i+1]<rint[i])&&(1.15*rint[i+1]<rint[i+2]))||((mod(rint[i]-rint[i+1])<0.3)&&(rint[i]<0.8&&rint[i+1]<0.8)&&(rint[i+2]>(1.2*avg(rint[i],rint[i+1]))))||(mod(rint[i+1]-rint[i+2])<0.3&&(rint[i+1]<0.8&&rint[i+2]<0.8)&&rint[i]>(1.2*avg(rint[i+1],rint[i+2])))){
			pvccount=1;
			vfcount=0;
			blockcount=0;
			start=i;
			i++;
			while(((1.15*rint[i+1]<rint[i])&&(1.15*rint[i+1]<rint[i+2]))||((mod(rint[i]-rint[i+1])<0.3)&&(rint[i]<0.8&&rint[i+1]<0.8)&&(rint[i+2]>(1.2*avg(rint[i],rint[i+1]))))||(mod(rint[i+1]-rint[i+2])<0.3&&(rint[i+1]<0.8&&rint[i+2]<0.8)&&rint[i]>(1.2*avg(rint[i+1],rint[i+2])))){
				pvccount++;
				i++;
				end=i;
				if(i>=llength-3){
					i=llength-3;
					end=i;
					break;
				}
			}
			if(pvccount>=3){
				cout<<"Ventricular tachycardia detected from "<<(double)rloc[start]/double(fs)<<" to "<<(double)rloc[end]/(double)fs<<endl;
				pvcflag=2;
			}
			else if(pvccount==2){
				cout<<"Ventricular couplet detected from "<<(double)rloc[start]/double(fs)<<" to "<<(double)rloc[end]/(double)fs<<endl;
				pvcflag=2;
			}
			else if(pvccount==1){
				pvcflag=1;
			}
		}


		////////////////////////////////////////////
		////////////Heart Block////////////////////
		//////////////////////////////////////////


		else if((2.2<rint[i+1]&&rint[i+1]<3.0)&&(mod(rint[i]-rint[i+1])<0.2||mod(rint[i+1]-rint[i+2])<0.2)){
			blockcount=1;
			vfcount=0;
			pvccount=0;
			pvcflag=2;
			start=i;
			i++;
			while((2.2<rint[i+1]&&rint[i+1]<3.0)&&(mod(rint[i]-rint[i+1])<0.2||mod(rint[i+1]-rint[i+2])<0.2)){
				blockcount++;
				end=i;
				if(i>=llength-3){
					i=llength-3;
					end=i;
					break;
				}
			}
			if(blockcount>1){
				cout<<"second degree heartblock detected from "<<(double)rloc[start]/double(fs)<<" to "<<(double)rloc[end]/(double)fs<<endl;
			}
		}
		else{
			pvcflag=0;
		}
		if(pos%2!=0&&pvcflag==1){
			if(pos==1){
				bistart=i;
			}
			bicount++;
			biendflag=0;
			pos++;
		}
		else if(pos%2==0&&pvcflag==0){
			bicount++;
			biendflag=0;
			pos++;
		}
		else{
			bicount=0;
			pos=1;
			biend=i;
			biendflag=1;
		}
		if(bicount>4&&biendflag==1){
			cout<<"bigeminy detected from "<<bistart<<" to "<<biend<<endl;
		}
	}
	system("pause");
	return 0;
}

void normalize(double arr[],int length){
	double max=0;
	int i=0;
	double val=0;
	max=((arr[0]<0)?-1*arr[0]:arr[0]);   ///assign modulus of ecg[0] to max
	for(i=1;i<length;i++){
		val=((arr[i]<0)?-1*arr[i]:arr[i]);          /////making negative values poistive (mod)
		if(val>max){
			max=val;
		}
	}
	for(i=0;i<length;i++){
		arr[i]/=max;
	}

}


void convolution(double arr[],double coef[],int length,int len2){
	int i=0,j=0;
	double sum=0;
	double temp[25000];
	for(i=0;i<(length+len2-1);i++){
		for(j=0;j<len2;j++){
			if((j-len2+1+i)>=length||(j-len2+1+i)<0){
				continue;
			}
			sum+=coef[j]*arr[j-len2+1+i];
	 	}
		temp[i]=sum;
			sum=0;
	}
	for(i=0;i<(length+len2-1);i++){
		arr[i]=temp[i];
	}
}

void square(double arr[],int length){
	int i=0;
	for(i=0;i<length;i++){
		arr[i]*=arr[i];
	}
}
double thresh(double arr[], int length){
	double max=0,sum=0,mean=0;
	int i=0;
	max=arr[0];   ///assign modulus of ecg[0] to max
	for(i=1;i<length;i++){
		if(arr[i]>max){
			max=arr[i];
		}
	}
	for(i=0;i<length;i++){
		sum+=arr[i];
	}
	mean=sum/(length);
	return (mean*max);
}

void interval(double ecg[],double thr,int lft[],int rgt [],int *llength,int *rlength,int length){
	int i=0;
	int lcount=0;
	int rcount=0;
	for(i=0;i<length;i++){
		if(ecg[i]>thr){
			
			lft[lcount]=i-38;
			lcount++;
			while(ecg[i]>thr){
				i++;
			}
			rgt[rcount]=i-38;
			rcount++;
		}
	}
	*llength=lcount;
	*rlength=rcount;
}
void maxima(double arr[], int left, int right, int* maxloc){
	int i=0;
	double max=arr[left];
	*maxloc=left;
	for(i=left+1;i<=right;i++){
		if(arr[i]>max){
			max=arr[i];
			*maxloc=i;
		}
	}
}
void minima(double arr[], int left, int right, int* minloc){
	int i=0;
	double min=arr[left];
	*minloc=left;
	for(i=left+1;i<=right;i++){
		if(arr[i]<min){
			min=arr[i];
			*minloc=i;
		}
	}
}
double pulserate(int length,int llength, int fs){
	//fs/=2;
	return ((double)llength*(double)fs*60)/(double)length;
}
void tachy_brady(int PR){
	if(100<PR&&PR<140){
		cout<<"sinus tachycardia"<<endl;
	}
	else if(140<PR&&PR<180){
		cout<<"Supra-ventricular tachycardia\n";
	}
	else if(PR>180){
		cout<<"ventricular tachycardia\n";
	}
	else if(60<PR&&PR<100){
		cout<<"heart beat normal\n";
	}
	else{
		cout<<"bradycardia\n";
	}
}
void bloodpressure(double PR,double* ss, double* ds){
	double map=0;
	map=(0.065*(double)PR*16.5)+5.0;
	*ds=map-10;
	*ss=(3*map)-(2*(*ds))+5;
	if(*ds<70||*ss<105){
		cout<<"Blood pressure is low\n";
	}
	else if(*ds>90||*ss>130){
		cout<<"Blood pressure is high\n";
	}
	else{
		cout<<"Blood pressure is normal\n";
	}
}
double avg(double a,double b){
	return ((a+b)/2.0);
}
double mod(double a){
	if(a>0){
		return a;
	}
	else{
		return (-1.0*a);
	}
}
