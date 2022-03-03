































default(parisize,"500M")
















\\define the parameters m=a*r
p_=17;
a_=9;
r_=5;








\\k is a finite field of cardinalitiy 17*9; f is a finite field of cardinalitiy 17*(9*5)
k=ffgen([p_,a_],'k);
f=ffgen([p_,a_*r_],'f);
m=ffembed(k,f);
kk=ffmap(m,k);












\\K_basis is one (absolute) basis w.r.t. F_17; F_basis_Relative is a relative basis of f w.r.t.k  
K_basis=vector(a_,i,kk^(i-1));
F_basis_Relative=vector(r_,i,f^(i-1));












\\define the Elliptic Curve y^2=x^3+1 on f; P_0 is a random point on E
E=ellinit([0,1],f);
\\P_0=random(E);
\\Q=ellgenerators(E);
P_0=[9*f^44 + 5*f^43 + 12*f^41 + 10*f^40 + 4*f^39 + 5*f^38 + 4*f^37 + 8*f^36 + 7*f^35 + 6*f^34 + 14*f^33 + 3*f^32 + 10*f^31 + 6*f^29 + 14*f^28 + 2*f^27 + 4*f^26 + 11*f^25 + 5*f^23 + 2*f^22 + 4*f^21 + 5*f^20 + 9*f^19 + 7*f^18 + 4*f^17 + 14*f^16 + 8*f^15 + 6*f^14 + 13*f^13 + 13*f^12 + 9*f^11 + 9*f^10 + 16*f^9 + 13*f^8 + 14*f^7 + 16*f^6 + 9*f^5 + 5*f^4 + 11*f^3 + 4*f + 3, 8*f^44 + 10*f^43 + 13*f^42 + 8*f^41 + 6*f^40 + 9*f^39 + 16*f^38 + 13*f^37 + 15*f^36 + 9*f^35 + 3*f^34 + f^33 + 8*f^32 + 13*f^31 + 4*f^30 + 14*f^29 + 5*f^28 + 8*f^27 + 12*f^26 + 13*f^25 + 11*f^24 + 11*f^23 + f^22 + 5*f^21 + 14*f^20 + 8*f^19 + 2*f^18 + 3*f^17 + 2*f^16 + 11*f^15 + 8*f^14 + f^13 + 16*f^12 + 3*f^11 + 11*f^10 + 6*f^9 + 9*f^8 + 8*f^7 + 9*f^6 + 6*f^5 + 10*f^4 + 14*f^3 + 16*f^2 + 13*f + 6];
Q=[[f^44 + 13*f^43 + 10*f^42 + 9*f^41 + 15*f^40 + 16*f^39 + 14*f^38 + 14*f^37 + 8*f^36 + f^35 + 14*f^34 + 14*f^33 + 6*f^32 + 13*f^31 + f^30 + 12*f^29 + 12*f^28 + 16*f^27 + 5*f^26 + 12*f^25 + 9*f^24 + 9*f^23 + 10*f^22 + 9*f^21 + 7*f^20 + 8*f^19 + 7*f^18 + 13*f^17 + 2*f^16 + 5*f^15 + 16*f^14 + 10*f^13 + 9*f^12 + 5*f^11 + 4*f^10 + 2*f^9 + 4*f^8 + 5*f^7 + 14*f^6 + 12*f^5 + 4*f^4 + 15*f^3 + 7*f^2 + 14*f + 4, 3*f^44 + 9*f^43 + 3*f^42 + 11*f^41 + 6*f^40 + 9*f^39 + 13*f^38 + 3*f^37 + 13*f^36 + 14*f^35 + 14*f^34 + 16*f^33 + 3*f^32 + 9*f^31 + 3*f^30 + 13*f^29 + 3*f^28 + 5*f^27 + 15*f^26 + f^25 + 9*f^24 + 12*f^23 + 10*f^22 + 8*f^21 + 4*f^20 + f^18 + 2*f^17 + 12*f^16 + 9*f^15 + 12*f^14 + 12*f^13 + 8*f^12 + 13*f^11 + 3*f^10 + 6*f^8 + 3*f^7 + 8*f^6 + 5*f^5 + 16*f^4 + 12*f^3 + 14*f + 8]];
















\\define the T operator
T(k_,P)=
{
     T_P=elladd(E,ellpow(E,P,k_),Q[1]);
     return (T_P);
}
































\\For z in f, ff in f, K_basis，define the composition of Phi and trace




phi_circ_trace(z,ff,K_basis,p,a)=
{
    summation=0;
    for(i=1,i=a,summation+=lift(trace(z*ff*K_basis[i]))/p^i);
    return(summation);
}
















\\For z in f, e_j go through F_basis_Relative, return the result of (Phi(Trace(z,e_j)))








Phi_Circ_Trace(z,F_basis_Rel,p,r,a)=
{
    res0=vector(r,i,phi_circ_trace(z,F_basis_Relative[i],K_basis,p,a));
    res0;
}
















R2R(qq)= qq+0.;
















F(kk,n)=
{  
     \\nn=fileopen("C:\\Users\\86159\\Desktop\\Res01\\Elliptic_random with n=" Str(n) ".txt", "a");
     RES=[];
     P=P_0;
     for(i=1,i=n,
        P=T(kk,P);
        Vcomp1=[R2R(x)|x<-(Phi_Circ_Trace(P[1],F_basis_Relative,p_,r_,a_))];
        Vcomp2=[R2R(x)|x<-(Phi_Circ_Trace(P[2],F_basis_Relative,p_,r_,a_))];
        res1=concat(Vcomp1,
                    Vcomp2);
        RES=concat(RES,[res1]);
                    \\filewrite(nn,Str(res1));
                    \\print(RES);
               
                    );
       
    \\fileclose(nn);
    return (RES);
}












FMatrix(kk,n)=
{  
     \\nn=fileopen("C:\\Users\\86159\\Desktop\\Res01\\Elliptic_random with n=" Str(n) ".txt", "a");
     RES=[];
     P=P_0;
     for(i=1,i=n,
        P=T(kk,P);
        Vcomp1=[R2R(x)|x<-(Phi_Circ_Trace(P[1],F_basis_Relative,p_,r_,a_))];
        Vcomp2=[R2R(x)|x<-(Phi_Circ_Trace(P[2],F_basis_Relative,p_,r_,a_))];
        res1=concat(Vcomp1,
                    Vcomp2);
        RES=concat(RES,[res1]~);
                    \\filewrite(nn,Str(res1));
                    \\print(RES);
               
                    );
       
    \\fileclose(nn);
    MRES=Mat(RES);
    return (MRES);
}




NEW_FMatrix(kk,LENGTH=10,HEIGHT)=
{  
     \\nn=fileopen("C:\\Users\\86159\\Desktop\\Res01\\Elliptic_random with n=" Str(n) ".txt", "a");
     RES2=[];
     P=P_0;
     for(i=1,i=HEIGHT,
        RES1=[];
        for(j=1,j=LENGTH,
            P=T(kk,P);
            Vcomp1=[R2R(x)|x<-(Phi_Circ_Trace(P[1],F_basis_Relative,p_,r_,a_))];
            Vcomp2=[R2R(x)|x<-(Phi_Circ_Trace(P[2],F_basis_Relative,p_,r_,a_))];
            VEC=concat(Vcomp1,
                     Vcomp2);
            RES1=concat(RES1,VEC);
            );
        MRes1=Mat(RES1);
        RES2=matconcat([Mat(RES2);MRes1]);
        );
       
    \\fileclose(nn);
    return (RES2);
}




\\-----------------------------------------------------------------------------------------




\\take the j-th col of M=[[,,],,[,,]]
Take_Data(M,j)=
{
   \\ for(i=1,i=length(M),print(M[i]));
    V=[];
    for(i=1,i=length(M),V=concat(V,M[i][j]));
    \\for(i=1,i=length(V),print(V[i]));
    return (V);
}




f_test(x)=
{
    if(x==0.0,return(0),
        return(1/x*cos(log(x)/x)));
}




Extract_positive(V)=
{
    R=[];
    for(i=1,i=length(V),if(V[i]>=0.0,R=concat(R,[V[i]])));
    return (R);
}
Extract_negative(V)=
{
    R=[];
    for(i=1,i=length(V),if(V[i]<0.0,R=concat(R,[V[i]])));
    return (R);
}
Integral_M(f,AA)=
{
    s=0;
    A=Take_Data(AA,2);
    fA=vector(length(A),i,f(A[i]));
    posi=Extract_positive(fA);
    \\print("Posi:    ",posi);
    nega=Extract_negative(fA);
    \\print("Nega:    ",nega);
    N=length(A);
    po_Sum=sum(i=1,i=length(posi),posi[i]);
    \\print("PoSum:   ", po_Sum);
    ne_Sum=sum(i=1,i=length(nega),nega[i]);
    s=po_Sum+ne_Sum;
    print(sqrt((s/N)^2));
    return(s/N);
}








F_test(V)=
{
    d=length(V);
    r=prod(i=1,i=d,exp(-V[i]/i));
    return (r);
}




Symbol_F_integral(d)=
{
    return (prod(i=1,d,i*(1-exp(-1/i))))
}
INTEGRAL_M100(F,M)=
{
   
    N=length(M);\\number of columns of M,
    print("N:=", N);
    mm=M[,^N];
    mt=mm[,^(N-1)];
    NN=length(M[,1]);
    S=sum(i=1,i=NN,F(mt[i,]));
    return (S/NN)
}




RINTEGRAL_M(F,M)=
{
   
    N=length(M);\\number of columns of M,
    print("N:=", N);
    \\mm=M[,^N];
    \\mt=mm[,^(N-1)];
    NN=length(M[,1]);
    S=sum(i=1,i=NN,F(M[i,]));
    return (S/NN)
}








\\Export the vector V,W to the txt;
Export_2_txt(V,W)=
{
    nn=fileopen("C:\\Users\\hrqin\\Desktop\\V,W " ".txt", "a");
    for(i=1,i=length(V),
        filewrite(nn,Str(V[i]) ", " Str(W[i]))
        );
    fileclose(nn);
    return(0);
}








EXPort_2_txt(M,i,j)=
{
    nn=fileopen("C:\\Users\\hrqin\\Desktop\\M" "[" Str(i) "]" "[" Str(j) "]" ".txt", "a");
    V1=Take_Data(M,i);
    print("V: ",V1);
    W1=Take_Data(M,j);
    print("W: ",W1);
    filewrite(nn,"V,W");
    for(k=1,k=length(V),
        print(V1[k],",   ",W1[k]);
        filewrite(nn,Str(V1[k]) ", " Str(W1[k]));
        );
    fileclose(nn);
    return(0);
}




MEXPort_2_txt(M,i,j)=
{
    nn=fileopen("C:\\Users\\hrqin\\Desktop\\Matrix" "[" Str(i) "]" "[" Str(j) "]" ".txt", "a");
    V1=M[,i];
    \\print("V: ",V1);
    W1=M[,j];
    \\print("W: ",W1);
    filewrite(nn,"V,W");
    for(k=1,k=length(V1),
        \\print(V1[k],",   ",W1[k]);
        filewrite(nn,Str(V1[k]) ", " Str(W1[k]));
        );
    fileclose(nn);
    return(0);
}




EXPORT_2_TXT(k1,Nu,i,j)=
{
    M=F(k1,Nu);
    nn=fileopen("C:\\Users\\hrqin\\Desktop\\Draw\\M" "[" Str(i) "]" "[" Str(j) "]" "(" Str(k1) "," Str(Nu) ")" ".txt", "a");
    V1=Take_Data(M,i);
   \\ print("V: ",V1);
    W1=Take_Data(M,j);
    \\print("W: ",W1);
    filewrite(nn,"V,W");
    for(k=1,k=length(V1),
        \\print(V1[k],",   ",W1[k]);
        filewrite(nn,Str(V1[k]) ", " Str(W1[k]));
        );
    fileclose(nn);
    return(0);
}








V_2_String(V)=
{
    string_sum="";
    for(i=1,i=length(V),
        strr=concat(Str(V[i])," ");
        string_sum=concat(string_sum,strr));
    return (string_sum);
}












new_MEXPORT_2_TXT(k1,Nu,L)=
{
    M=NEW_FMatrix(k1,L/10,Nu);
    nn=fileopen("C:\\Users\\Admin\\Documents\\M"  "(k=" Str(k1) ", N=" Str(Nu) ", L=" Str(L) ")" ".txt", "a");




    for(k=1,k=Nu,
        filewrite(nn,V_2_String(M[k,]));
        );
    fileclose(nn);
    return(0);
}


Direct_MEXPORT_2_TXT(k1,Nu,L)=
{
    M=NEW_FMatrix(k1,L/10,Nu);
    nn=fileopen("C:\\Users\\Admin\\Documents\\M"  "(k=" Str(k1) ", N=" Str(Nu) ", L=" Str(L) ")" ".txt", "a");
    for(i=1,i=Nu,
        for(j=1,j=L,
            filewrite(nn,M[i,j] )
        )
    );
    fileclose(nn);
    return(0);
}




FM_export(kk,LENGTH=10,HEIGHT)=
{  
     nn=fileopen("C:\\Users\\Admin\\Documents\\M"  "(k=" Str(kk) ", H=" Str(HEIGHT) ", L=" Str(LENGTH) ")" ".txt", "a");
     \\localprec(4);
     P=P_0;
     for(i=1,i=HEIGHT,
        RES1=[];
        for(j=1,j=LENGTH,
            P=T(kk,P);
            Vcomp1=[R2R(x)|x<-(Phi_Circ_Trace(P[1],F_basis_Relative,p_,r_,a_))];
            Vcomp2=[R2R(x)|x<-(Phi_Circ_Trace(P[2],F_basis_Relative,p_,r_,a_))];
            VEC=concat(Vcomp1,
                     Vcomp2);
            RES1=concat(RES1,VEC);
            );
        for(j=1,j=length(RES1),
            filewrite(nn,RES1[j])
            );
        fileflush(nn);
        );
       
    fileclose(nn);
   
    \\localprec(default(realprecision));
}


fm_export(kk,LENGTH=10,HEIGHT)=
{  
     default(logfile,"C:\\Users\\Admin\\Documents\\M"  "(k=" Str(kk) ", H=" Str(HEIGHT) ", L=" Str(LENGTH) ")" ".log");
     default(log,1);
     P=P_0;
     for(i=1,i=HEIGHT,
        RES1=[];
        for(j=1,j=LENGTH,
            P=T(kk,P);
            Vcomp1=[R2R(x)|x<-(Phi_Circ_Trace(P[1],F_basis_Relative,p_,r_,a_))];
            Vcomp2=[R2R(x)|x<-(Phi_Circ_Trace(P[2],F_basis_Relative,p_,r_,a_))];
            VEC=concat(Vcomp1,
                     Vcomp2);
            RES1=concat(RES1,VEC);
            );
        for(j=1,j=length(RES1),
            if((i<HEIGHT) || (j<length(RES1)),
                printf("%1.10f\n",RES1[j]),
                printf("%1.10f",RES1[j]))
            );
        );


}
\\-----------------------------------------------------------------------------------------
aff(x,L0)=
{
    return (L0*x);
}




AFF(M,ind1,ind2,L1,L2)=
{
    nn=fileopen("C:\\Users\\hrqin\\Desktop\\Draw\\aff_M" "[" Str(ind1) "]" "[" Str(ind2) "]" "(L1=" Str(L1) ", L2=" Str(L2) ")" ".txt", "a");
    M1=Take_Data(M,ind1);
    M2=Take_Data(M,ind2);




    V1=vector(i=length(M1),i,aff(M1[i],L1));
    V2=vector(i=length(M2),i,aff(M2[i],L2));
    print(V1);
    print(V2);




    filewrite(nn,"V,W");
    for(k=1,k=length(V1),
        filewrite(nn,Str(V1[k]) ", " Str(V2[k]));
        );
    fileclose(nn);
}




EX_AFF(M,ind1,ind2,x1,x2,y1,y2)=
{
    nn=fileopen("C:\\Users\\hrqin\\Desktop\\Draw\\aff_M" "[" Str(ind1) "]" "[" Str(ind2) "]" "(L1=" Str(x2-x1) ", L2=" Str(y2-y1) ")" ".txt", "a");
    L1=x2-x1;
    L2=y2-y1;
    M1=Take_Data(M,ind1);
    M2=Take_Data(M,ind2);




    V1=vector(i=length(M1),i,aff(M1[i],L1)+x1);
    V2=vector(i=length(M2),i,aff(M2[i],L2)+y1);
    print(V1);
    print(V2);




    filewrite(nn,"V,W");
    for(k=1,k=length(V1),
        filewrite(nn,Str(V1[k]) ", " Str(V2[k]));
        );
    fileclose(nn);
}
\\-----------------------------------------------------------------------------------------












Estimate_Bound_func_1_dim(f,x1,x2,Density=1000)=
{
    \\Generating the uniformly distributed sequence of points.
    \\Density=1000;
    M=F(1,Density);
    M1=Take_Data(M,1);
    V0=vector(i=Density,i,aff(M1[i],x2-x1)+x1);




    VMax=vecmax([f(x)|x<-V0]);
    MaxV=vecmax([f(x1),f(x2),VMax]);
    \\print([f(x1),f(x2),VMax]);
    CMV=R2R(ceil(MaxV)) ;
    return (CMV);
}
f_sin_test(x)=
{
    if(x==0,0,sin(1/x)+1)
}
Monte_Carlo_Ingeral_Positive_2D(func,TN,x1,x2)=
{
    C=0;
    YMax=Estimate_Bound_func_1_dim(func,x1,x2);
    M=F(1,TN);
    X0=Take_Data(M,1);
    X=vector(i=length(X0),i,aff(X0[i],x2-x1)+x1);
    Y0=Take_Data(M,2);
    Y=vector(i=length(Y0),i,aff(Y0[i],YMax));
    for(i=1,i=TN,
        if(func(X[i])>=Y[i],C++);
    );
    return (R2R(C/TN)*(YMax*(x2-x1)));
}
MC_Integral_Positive_2D_Elementary(func,TN)=
{
    C=0;
    M=F(1,TN);
    X0=Take_Data(M,1);
    X=vector(i=length(X0),i,aff(X0[i],1));
    Y0=Take_Data(M,2);
    Y=vector(i=length(Y0),i,aff(Y0[i],2));
    for(i=1,i=TN,
        if(func(X[i])>=Y[i],C++);
    );
    return (R2R(C/TN)*(2*1));
}




N_plus_1_MC(func,M)=
{
    C=0;
    N=length(M);
    mm=M[,^N];
    print("N=",N);
    X0=mm[,^(N-1)];
    Y0=M[,(N-1)];
    NN=length(M[,1]);
    print("NN=",NN);
    for(i=1,i=NN,
        if(func(X0[i,])>=Y0[i],C++);
        );
    return (R2R(C/NN));
}
\\M=F(5000);
\\M1=Take_Data(M,1);
\\M2=Take_Data(M,2);
\\s=plothsizes();
\\plotinit(0,s[2]-1,s[2]-1);
\\plotscale(0,-1,1,-1,1);
\\plotpoints(0,Take_Data(M,1),Take_Data(M,2));
\\plotdraw(0);






























