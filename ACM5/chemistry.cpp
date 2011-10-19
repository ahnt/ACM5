#include "chemistry.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>


t_chemistry::t_chemistry()
{

}

t_chemistry::~t_chemistry()
{
}

void t_chemistry::setup(void)
{
	int32_t i;
	unsigned char k;
	int32_t z=0;
	int32_t x,y;
	int u=0;
	unsigned char ca,cb;
	int da,db,dwo_a,dwo_b,dap,dbp;
	FILE *f;
	printf("	* INITIALIZING CHEMISTRY...");
	for(i=0;i<0x1000000;i++)
	{
		if(this->is_valid_compound(i)==0)
		{
			this->compound[z]=i;
			for(k=0;k<4;k++)
				this->compound_content[z][k]=0;
			for(k=0;k<12;k++)
			{
				this->compound_content[z][(i>>(k*2))&3]++;
			}
			z++;
		}
	}
	/*
	f=fopen("c:\\aff3.aff","r+b");
	for(z=0;z<0x1000000;z++)
	{
		fscanf(f,"%c%c",&ca,&cb);
		q_aff[z&0xffffff]=ca+(cb*256);
	}
	fclose(f);
	*/
	printf("	...COMPLETED *\n");
	/*
		for(x=0;x<mx;x++)
			for(y=0;y<my;y++)
			{
				for(z=0;z<608;z++)
				{
					this->n[x][y][z]=(double) start_compounds;
				}
				n[x][y][608]=(double) start_compounds*608;
			}*/
	//this->np[608]=(double) start_compounds*(double) 608;
	printf("	* LOADING POSSIBLE CHEMICAL REACTIONS...");
	this->q_pcr=new t_pcr[max_pcr];
	f=fopen("pos_chem_reactions.txt","r+t");
	for(z=0;z<max_pcr;z++)
	{
		fscanf(f,"%i	%i	%i	%i\n",&da,&db,&dwo_a,&dwo_b);
		this->q_pcr[z].A=(int16_t) da;
		this->q_pcr[z].B=(int16_t) db;
		this->q_pcr[z].wo_A=(unsigned char)(dwo_a&15);
		this->q_pcr[z].wo_B=(unsigned char)(dwo_b&15);
		if((da<105)&&(db<105)) u++;
		/*
		if((da==6)&&(db==11))
		{
			printf("\n");
			printf("%i	%i	%i	%i	%i\n",z,da,db,dwo_a,dwo_b);
			this->quick_split(da,db,dwo_a,dwo_b,dap,dbp);
			this->show_str_compound(compound[da]);
			printf("	");
			this->show_str_compound(compound[db]);
			printf("	");
			this->show_str_compound(dap);
			printf("	");
			this->show_str_compound(dbp);
			printf("\n");
		}
		*/
	}
	fclose(f);
	printf("	...DONE *\n");
	printf("	%i reaktions below 54\n",u);
}
double t_chemistry::affinity(int32_t a, int32_t b)
{
	double r=0.0;
	double d;
	int z;
	for(z=0;z<12;z++)
	{
		d=(((a>>(z*2))&3)^((b>>(z*2))&3));
		r=r+(d*d);
	}
	return (double)1.0-sqrt((double)r/(double)192);
}
int32_t t_chemistry::is_valid_compound(int32_t to_test)
{
	int32_t i=0;
	int32_t t=0;
	int32_t lc=0;
	if((to_test&3)!=0)
	{
		while(t<12)
		{
			if((t&1)==0) i+=(to_test>>(t*2))&3;
			else i-=(to_test>>(t*2))&3;
			if(t<(11)) 
			{
				lc=t&1;
				if((i==0)&&(((to_test>>((t+1)*2))&3)!=0)) i=1000;
				if((i!=0)&&(((to_test>>((t+1)*2))&3)==0)) i=1000;
				if(((t&1)==0)&&(i<0)) i=1000;
				if(((t&1)==1)&&(i>0)) i=1000;
			}
			t++;
		}
	}
	else
	i=1;
	if(to_test<0) i=1;
	return i;
}
unsigned char t_chemistry::contains_site(int32_t A, int32_t B)
{
	unsigned char mA[12],mB[12];
	unsigned char x,y,s,d;
	unsigned char z;
	for(x=0;x<12;x++)
	{
		mA[x]=(A>>(x*2))&3;
		mB[x]=(B>>(x*2))&3;
	}
	z=0;
	s=0;
	for(x=0;x<12;x++)
	{
		d=0;
		for(y=0;y<(12-x);y++)
			if(mA[x+y]==mB[y]) d++;
		if (d>s)
		{
			s=d;
			z=x;
		}
	}
	if(s!=0)
	return z;
	else 
	return -1;
}

unsigned char t_chemistry::split(int32_t O1,int32_t O2,int32_t O3,int32_t O4,int32_t &O5,int32_t &O6)
{
	unsigned char mO1[12],mO2[12],mO3[12],mO4[12];
	int16_t x,y,z,s,d;
	int sa,sb;
	int i=0;
	for(x=0;x<12;x++)
	{
		mO1[x]=(O1>>(x*2))&3;
		mO2[x]=(O2>>(x*2))&3;
		mO3[x]=(O3>>(x*2))&3;
		mO4[x]=(O4>>(x*2))&3;
	}
	s=0;
	sa=0;
	for(x=0;x<12;x++)
	{
		d=0;
		for(y=0;y<(12-x);y++)
			if(mO1[x+y]==mO3[y]) d++;
		if (d>s)
		{
			s=d;
			sa=x;
		}
	}
	s=0;
	sb=0;
	for(x=0;x<12;x++)
	{
		d=0;
		for(y=0;y<(12-x);y++)
			if(mO2[x+y]==mO4[y]) d++;
		if (d>s)
		{
			s=d;
			sb=x;
		}
	}
	O5=0;
	O6=0;
	for(z=11;z>-1;z--)
	{
		if(z<sa) O5=(O5<<2)+mO1[z]; else O5=(O5<<2)+mO2[sb+(z-sa)];
		if(z<sb) O6=(O6<<2)+mO2[z]; else O6=(O6<<2)+mO1[sa+(z-sb)];
	}
	if ((is_valid_compound(O5)==0)&&(is_valid_compound(O6)==0))
	return 1;
	else 
	return 0;
}

void t_chemistry::show_str_compound(int32_t i)
{
	int z;
	for(z=0;z<12;z++)
		printf("%i",(i>>(z*2))&3);
}

unsigned char t_chemistry::length(int32_t O)
{
	unsigned char i=0;
	int z;
	for(z=0;z<12;z++)
		if(((O>>(z*2))&3)!=0) i++;
	return i;
}

unsigned char t_chemistry::delta_e(int32_t O1,int32_t O2,int32_t O3,int32_t O4)
{
	return abs(length(O1)-length(O2))-abs(length(O3)-length(O4));
}

double t_chemistry::c(double x, double y, int32_t comp)
{
	double d;
	d=sqrt(((x-xp[comp])*(x-xp[comp]))+((y-yp[comp])*(y-yp[comp])))/50;
	d=powf(M_E,(-1*d*d)/2)/sqrt(2*M_PI)*(np[comp]/np[608]);
	if(sqrt((x*x)+(y*y))>(double)200)
	{
		d=(double)0.0;
	}
	return d;
}

void t_chemistry::show_concentrations(void)
{
	int z;
	for(z=0;z<608;z++)
	{
		printf("%i	",z);
		this->show_str_compound(z);
		printf("	%f\n",np[z]);
	}
}

unsigned char t_chemistry::quick_split(int A, int B, int wo_A, int wo_B, int &AP, int &BP)
{
	unsigned char mO1[12],mO2[12],mO3[24],mO4[24];
	int16_t x,y,z,s,d;
	int sa,sb;
	int i=0;
	for(x=0;x<12;x++)
	{
		mO1[x]=(this->compound[A]>>(x*2))&3;
		mO2[x]=(this->compound[B]>>(x*2))&3;
	}
	AP=0;
	BP=0;
	for(x=0;x<12;x++)
	{
		if(x<wo_A) mO3[x]=mO1[x];
		else
		{
			if((wo_B+x-wo_A)<12) mO3[x]=mO2[wo_B+x-wo_A];
			else mO3[x]=0;
		}
		if(x<wo_B) mO4[x]=mO2[x];
		else
		{
			if ((wo_A+x-wo_B)<12) mO4[x]=mO1[wo_A+x-wo_B];
			else mO4[x]=0;
		}
	}
	/*
	for(x=0;x<12;x++)
		printf("%i",mO1[x]);
	printf("\n");
	for(x=0;x<12;x++)
		printf("%i",mO2[x]);
	printf("\n");
	for(x=0;x<12;x++)
		printf("%i",mO3[x]);
	printf("\n");
	for(x=0;x<12;x++)
		printf("%i",mO4[x]);
	printf("\n");
	printf("\n");*/
	for(x=11;x>-1;x--)
	{
		AP=(AP<<2)+(mO3[x]);
		BP=(BP<<2)+(mO4[x]);
	}
	if ((is_valid_compound(AP)==0)&&(is_valid_compound(BP)==0))
	return 1;
	else 
	return 0;
}

int t_chemistry::inttolinked(int v)
{
	int z=0;
	while(this->compound[z]!=v) z++;
	if(z>608) printf("PANIK\n");
	return z;
}