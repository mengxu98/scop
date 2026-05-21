// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X){
  int i, p=X.cols(), N=X.rows();
  Eigen::VectorXd mX(p), sdX(p);
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  return List::create(Named("x")=X, Named("sd")=sdX, Named("m")=mX);
}

/*****  Omega  *****/
// [[Rcpp::export]]
List OmegaC(Eigen::MatrixXd & Omega, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);

  //Omega.diagonal().setZero();
  Eigen::SparseMatrix<double> OmegaS=Omega.sparseView();

  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }

  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }

  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}

/*****  Sparse Omega  *****/
// [[Rcpp::export]]
List OmegaSC(Eigen::SparseMatrix<double> & OmegaS, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);

  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }

  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }

  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}





/////////////////////////////////
/////   Linear Regression   /////
/////////////////////////////////

/*****  LM: Lambda path (max) inner product <xj,y> *****/
// [[Rcpp::export]]
double maxLambdaLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd wbeta, int N0, int p){
  int i;
  double LiMax=0.0, LiMaxi=0.0;

  for (i=0; i<p; ++i) {
    if (wbeta(i) > 0.0) {
      LiMaxi=std::abs(y.transpose()*X.col(i))/wbeta(i); // <xj,y>/N0
      if (LiMaxi > LiMax) {
        LiMax=LiMaxi;
      }
    }
  }

  LiMax=LiMax/N0/alpha;

  return(LiMax);
}



/*****  Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvTrimLmC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco,
                          Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double a0) {
  int i, j;
  Eigen::VectorXd RSS, xbF=Eigen::VectorXd::Zero(NF); // xb=Eigen::VectorXd::Zero(N),

  yF=yF.array()-a0;

  if(nn2>0){
    RSS.setZero(nn2); //nn2= # of part of data

      if(nn==0){
        RSS(0)=yF.squaredNorm();
      }else{
        for(i=0;i<nn;i++){
          j=loco(i); //   index of nonzero beta
          xbF+=XF.col(j)*beta(i); //
            RSS(i)=(yF-xbF).squaredNorm();
        }
      }

    if(nn2>nn && nn>0){
      for(i=nn;i<nn2;i++){RSS(i)=RSS(nn-1);}
    }

    if(nn2>nn && nn==0){
      for(i=nn;i<nn2;i++){RSS(i)=RSS(0);}
    }

  }else{
    RSS.setZero(1);
    RSS(0)=yF.squaredNorm();
  }

  return(RSS);
}


/*****  LM: Enet (L1+L2)  *****/
// [[Rcpp::export]]
List EnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
             double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
             int p, int N0, double thresh, int maxit, double thresh2){

  int  i, j, it=0, il, iadd, ia=0;
  double lambda2, zi, obj0, obj1, b0, db0, objQi, objQj, rss0, rss1, RSS0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd RSS=Eigen::VectorXd::Zero(nlambda), RSQ=Eigen::VectorXd::Zero(nlambda);
  double xr, dbMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }

  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          if(zi>lambda1(j)){
            // b0=(zi-lambda1(j))/(lambda2*wbeta(j)+1); // x*x/N=1
            b0=(zi-lambda1(j))/(lambda2+1.0); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1(j)){
            // b0=(zi+lambda1(j))/(lambda2*wbeta(j)+1);
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N0/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}


    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("rsq")=RSQ,
                      Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}


/*****  LM: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
               double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
               int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double thresh2){

  int i, j, it=0, il, iadd, ia=0;
  double lambda2, zi, obj0, obj1, rss0, rss1, b0, db0, objQi, objQj, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda), BetaSTD=Eigen::MatrixXd::Zero(p,nlambda); // beta matrix for different lambdas
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXd RSSp(nlambda), RSS(nlambda), RSQ(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  double xr, dbMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);
  double a0=0.0, my=0.0;
  Eigen::MatrixXd predY=Eigen::MatrixXd::Zero(NF, nlambda);

  Eigen::VectorXd a0S=Eigen::VectorXd::Zero(nlambda);


  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
//    X.col(i)/=sdX(i);
  }
  my=y.mean();
  y=y.array()-my;

  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1.0); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    a0=my; xbF.setZero(NF);
    for(i=0;i<ia;i++){
      j=active(i);
      xbF+=XF.col(j)*Beta(j,il);
      a0-=mX(j)*Beta(j,il);
    }
    xbF=xbF.array()+a0;

    predY.col(il)=xbF;
    RSSp(il)=(yF-xbF).squaredNorm();
    a0S(il)=a0;

    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("predY")=predY, Named("a0S")=a0S,
                      Named("RSS")=RSS, Named("rsq")=RSQ, Named("RSSp")=RSSp, Named("nlambda")=il));
}



/*****  LM: Network (L1+La)  *****/
  // [[Rcpp::export]]
List NetLmC(Eigen::MatrixXd & X, Eigen::VectorXd & y, double alpha,
            Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
            Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
            int p, int N0, double thresh, int maxit, double thresh2){

  int i, j, ij, m, it=0, il, iadd, ia=0;
  double lambda2, zi, zi2, objQi=0.0, objQj, obj0, obj1, rss0, rss1, b0, db0, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N0);
  Eigen::VectorXd RSS=Eigen::VectorXd::Zero(nlambda), RSQ(nlambda);
  double xr, dbMax, lambdaMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }

  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);} // Omega: w_kl/sqrt(d_k*d_l),L=I-Omega ; L=SLS (included sign of beta)
          }
          zi+=lambda2*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2); //  beta^T*L*beta
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N0/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss1-rss0)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){
          flag(il)=2;
          goto exit;}
        if(it>=maxit){
          flag(il)=1;
          break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i);
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)=std::abs(y.dot(X.col(i))/N0+lambda2*zi2);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
    return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("rsq")=RSQ,
                        Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  LM: Network (L1+La) cross-validation *****/
// [[Rcpp::export]]
List cvNetLmC(Eigen::MatrixXd & X, Eigen::VectorXd & y,double alpha,
              Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
              Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
              int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double thresh2){

  int i, j, ij, m, it=0, il, iadd, ia=0;
  double lambda2, zi, zi2, objQi=0.0, objQj, obj0, obj1, rss0, rss1, b0, db0, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  Eigen::VectorXd RSS(nlambda), RSQ(nlambda), RSSp(nlambda);
  double xr, dbMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  double a0=0.0, my=0.0;
  Eigen::MatrixXd predY=Eigen::MatrixXd::Zero(NF, nlambda);

  Eigen::VectorXd a0S=Eigen::VectorXd::Zero(nlambda);

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
//    X.col(i)/=sdX(i);
  }
  my=y.mean();
  y=y.array()-my;

  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);} // Omega: w_kl/sqrt(d_k*d_l),L=I-Omega ; L=SLS (included sign of beta)
          }
          zi+=lambda2*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2); //  beta^T*L*beta
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss1-rss0)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i);
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)=std::abs(y.dot(X.col(i))/N+lambda2*zi2);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    a0=my; xbF.setZero(NF);
    for(i=0;i<ia;i++){
      j=active(i);
      xbF+=XF.col(j)*Beta(j,il);
      a0-=mX(j)*Beta(j,il);
    }
    xbF=xbF.array()+a0;
    predY.col(il)=xbF;

    RSSp(il)=(yF-xbF).squaredNorm();
    a0S(il)=a0;

    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("predY")=predY, Named("a0S")=a0S,
                      Named("RSS")=RSS, Named("RSSp")=RSSp, Named("rsq")=RSQ, Named("nlambda")=il));
}





///////////////////
/////   Cox   /////
///////////////////

/*****  Cox: Lambda path (max)  *****/
// [[Rcpp::export]]
double maxLambdaCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, int N,
                     Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
                     int n, double alpha, Eigen::VectorXd wbeta, int N0, int p){
  int i, j, q;
  double denS=N, c1=0.0;
  Eigen::VectorXd lli(N);
  double LiMax=0.0, LiMaxi=0.0;

  for(i=0;i<n;i++){
    c1+=(nevent1(i)/denS);denS-=nevent(i);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){lli(j)=tevent(j)-c1;}
  }


  for (i=0; i<p; ++i) {
    if (wbeta(i) > 0.0) {
      LiMaxi=std::abs(lli.transpose()*X.col(i))/wbeta(i); // <xj,y>/N0
      if (LiMaxi > LiMax) {
        LiMax=LiMaxi;
      }
    }
  }

  LiMax=LiMax/N0/alpha;

  return(LiMax);
}



/*****  Derivatives of log-pl of eta (1st&2nd order),  ties  *****/
void dletaCm(Eigen::VectorXd& exb, Eigen::VectorXd& tevent, int& N,
               Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1,
               int& n, Eigen::VectorXd& pl1, Eigen::VectorXd& pl2, int& ifast, int& itwo){
    int i, j, q, ipl2=0;
    double denSi, c1=0.0, c2=0.0;
    Eigen::VectorXd denS(n);

    if(ifast==0 || itwo==1)goto two;
    denSi=exb.sum();
    for(i=0;i<n;++i){
      c1+=(nevent1(i)/denSi);c2+=(nevent1(i)/pow(denSi, 2));
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
        denSi-=exb(j);
        pl1(j)=tevent(j)-exb(j)*c1;
        pl2(j)=exb(j)*(c1-exb(j)*c2);
        if(pl2(j)<=0.0)ipl2=1;
      }
    }
    if(ipl2==1){itwo=1;if(ifast==0){goto two;}}
    return;

    two:
      denSi=0.0;c1=0.0;c2=0.0;
    for(i=n-1;i>=0;--i){
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){denSi+=exb(j);}
      denS(i)=denSi;
    }
    for(i=0;i<n;++i){
      c1+=(nevent1(i)/denS(i));c2+=(nevent1(i)/pow(denS(i), 2));
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
        pl1(j)=tevent(j)-exb(j)*c1;
        pl2(j)=exb(j)*(c1-exb(j)*c2);
    }
  }
}

/*****  Log-pl of eta,  ties  *****/
// [[Rcpp::export]]
double pletaCm(Eigen::VectorXd& xb, Eigen::VectorXd& exb, Eigen::VectorXi& nevent,
               Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n, int& ifast, int& itwo){
  int i, j, q, iSS=0;
  double ll=0.0, SSi;
  Eigen::VectorXd SS(n);

  if(ifast==0 || itwo==1)goto two;
  SSi=exb.sum();
  for(i=0;i<n;++i){
    if(SSi<=0.0)iSS=1;
    for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){ll+=xb(j);}
    ll-=nevent1(i)*log(SSi);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi-=exb(j);}
  }
  if(iSS==1){itwo=1;if(ifast==0){goto two;}}
  return(ll);

  two:
    ll=0.0;SSi=0.0;
    for(i=n-1;i>=0;--i){
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi+=exb(j);}
      SS(i)=SSi;
    }
    for(i=0;i<n;++i){
      for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){
        ll+=xb(j)-log(SS(i));
      }
    }
    return(ll);
}


/*****  Cox: Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvTrimCoxC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco,
                           Eigen::MatrixXd XF, int NF,
                           Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF,
                           Eigen::MatrixXd X, int N,
                           Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n,
                           int ifast, int itwo){
  int i, j;
  double lli, lfi;
  Eigen::VectorXd cv, xb=Eigen::VectorXd::Zero(N), xbF=Eigen::VectorXd::Zero(NF);
  Eigen::VectorXd exb(N), exbF(NF);

  if(nn2>0){
    cv.setZero(nn2);

    if(nn==0){
      exb=(xb.array()).exp();
      lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      exbF=(xbF.array()).exp();
      lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
      cv(0)=lfi-lli;

    }else{
      for(i=0;i<nn;i++){
        j=loco(i);
        xb+=X.col(j)*beta(i);exb=(xb.array()).exp();
        lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        xbF+=XF.col(j)*beta(i);exbF=(xbF.array()).exp();
        lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
        cv(i)=lfi-lli;
      }
    }

    if(nn2>nn && nn>0){
      for(i=nn;i<nn2;i++){cv(i)=cv(nn-1);}
    }

    if(nn2>nn && nn==0){
      for(i=nn;i<nn2;i++){cv(i)=cv(0);}
    }

  }else{
    cv.setZero(1);

    exb=(xb.array()).exp();
    lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
    exbF=(xbF.array()).exp();
    lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
    cv(0)=lfi-lli;
  }

  return(cv);
}


/*****  Cox: Enet (L1+L2)  *****/
  // [[Rcpp::export]]
List EnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent,
              double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
              int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
              int n, int p, int N0, double thresh, int maxit, int ifast){

  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda),BetasSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double lambdaMax;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0;objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();
  }//for lambda

  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("lambda")=lambda, Named("ll")=locbeta, Named("nlambda")=il));
}


/*****  Cox: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent,
                double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
                int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
                int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF,
                int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){

  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda),BetasSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double mxi;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);

    mxi=XF.col(i).mean();
    XF.col(i)=XF.col(i).array()-mxi;
  }

  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0;objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();

    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*Betas(j,il);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda

  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}



/*****  Cox: Network (L1+La)  *****/
  // [[Rcpp::export]]
List NetCoxC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha,
             Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
             Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
             int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
             int n, int p, int N0, double thresh, int maxit, int ifast){

  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda), BetasSTD=Eigen::MatrixXd::Zero(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double lambdaMax;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2i*zi2;

          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();
  }//for lambda

  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("ll")=locbeta, Named("nlambda")=il));
}


/*****  Cox: Network (L1+La)  cross-validation  *****/
// [[Rcpp::export]]
List cvNetCoxC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha,
               Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
               Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
               int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
               int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF,
               int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){

  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda), BetasSTD=Eigen::MatrixXd::Zero(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double mxi;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);

    mxi=XF.col(i).mean();
    XF.col(i)=XF.col(i).array()-mxi;
  }

  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2i*zi2;

          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();

    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*Betas(j,il);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda

  exit:
    if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}





///////////////////////////////////
/////   Logistic Regression   /////
///////////////////////////////////

/*****  Log: Lambda path (max) inner product <xj,y> *****/
// [[Rcpp::export]]
double maxLambdaLogC(Eigen::MatrixXd X, Eigen::VectorXd Z,
                     double alpha, Eigen::VectorXd wbeta, int N0, int p){
  int i;
  double LiMax=0.0, LiMaxi=0.0;

  for (i=0; i<p; ++i) {
    if (wbeta(i) > 0.0) {
      LiMaxi=std::abs(Z.transpose()*X.col(i))/wbeta(i); // <xj,y>/N0
      if (LiMaxi > LiMax) {
        LiMax=LiMaxi;
      }
    }
  }

  LiMax=LiMax/N0/alpha;

  return(LiMax);
}


/*****  Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvTrimLogC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco,
                           Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double threshP) {
  int i, j;
  Eigen::VectorXd Dev(nn2), xbF=Eigen::VectorXd::Zero(NF); // xb=Eigen::VectorXd::Zero(N),
  Eigen::ArrayXd p0(NF);

  for(i=0;i<nn;i++){
    j=loco(i); //   index of nonzero beta
    xbF+=XF.col(j)*beta(i);

    for (j=0; j<NF; ++j) {
      p0(j)=1.0/(1.0+exp(-xbF(j)));
      if (p0(j) < threshP) {
        p0(j)=threshP;
      } else if (p0(j) > (1.0-threshP)) {
        p0(j)=1.0-threshP;
      }
    }

    Dev(i)=(yF.array()*p0.log()+(1.0-yF.array())*(1.0-p0).log()).sum()*(-2.0);
  }


  if(nn2>nn && nn>0){
    for(i=nn;i<nn2;i++){Dev(i)=Dev(nn-1);}
  }

  return(Dev);
}


/*****  Log: Enet (L1+L2)  *****/
// [[Rcpp::export]]
List EnetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
              double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
              int p, int N0, double thresh, int maxit, double threshP){

  int  i, j, i2, it=0, il, iadd, ia=0;
  double zi, b0, db0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;


  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  // Lambda path
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLogC(X, Z, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLogC(X, Z, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }


          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update


        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

    iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(Z.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));

  }//for lambda

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  Log: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
                double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
                int p, int N0, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double threshP){

  int  i, j, i2, it=0, il, iadd, ia=0;
  double zi, b0, db0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;

  Eigen::ArrayXd xbF(NF), pF0(NF), LLF(nlambda);
  double llF0;

  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }

          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

      iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(Z.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));


    // Predict Deviance
    xbF=XF*Beta.col(il);
    for (i2=0; i2<NF; ++i2) {
      pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
      if (pF0(i2) < threshP) {
        pF0(i2)=threshP;
      } else if (pF0(i2) > (1.0-threshP)) {
        pF0(i2)=1.0-threshP;
      }
    }
    LLF(il)=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  }//for lambda

  xbF=log(yF.mean()/(1.0-yF.mean()))*XF.col(0).array();
  for (i2=0; i2<NF; ++i2) {
    pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
    if (pF0(i2) < threshP) {
      pF0(i2)=threshP;
    } else if (pF0(i2) > (1.0-threshP)) {
      pF0(i2)=1.0-threshP;
    }
  }
  llF0=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("LLF")=LLF.head(il), Named("llF0")=llF0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}





/*****  Log: Network (L1+La)  *****/
// [[Rcpp::export]]
List NetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
             double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
             Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
             int p, int N0, double thresh, int maxit, double threshP){

  int  i, j, i2, m, ij, it=0, il, iadd, ia=0;
  double zi, zi2, b0, db0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;


  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  // Lambda path
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLogC(X, Z, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLogC(X, Z, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          zi2=0.0;
          for(ij=0; ij<nadj(j); ++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2(j)*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }

          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

    iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){

        di(i)=std::abs(Z.dot(X.col(i))/N0);

        zi2=0.0;
        for(ij=0; ij<nadj(i); ++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)+=lambda2(i)*zi2;

        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));

  }//for lambda

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  Log: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvNetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
               double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
               Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
               int p, int N0, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double threshP){

  int  i, j, i2, ij, m, it=0, il, iadd, ia=0;
  double zi, zi2, b0, db0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;

  Eigen::ArrayXd xbF(NF), pF0(NF), LLF(nlambda);
  double llF0;

  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          zi2=0.0;
          for(ij=0; ij<nadj(j); ++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2(j)*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }

          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

      iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){

        di(i)=std::abs(Z.dot(X.col(i))/N0);

        zi2=0.0;
        for(ij=0; ij<nadj(i); ++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)+=lambda2(i)*zi2;

        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));


    // Predict Deviance
    xbF=XF*Beta.col(il);
    for (i2=0; i2<NF; ++i2) {
      pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
      if (pF0(i2) < threshP) {
        pF0(i2)=threshP;
      } else if (pF0(i2) > (1.0-threshP)) {
        pF0(i2)=1.0-threshP;
      }
    }
    LLF(il)=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  }//for lambda

  xbF=log(yF.mean()/(1.0-yF.mean()))*XF.col(0).array();
  for (i2=0; i2<NF; ++i2) {
    pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
    if (pF0(i2) < threshP) {
      pF0(i2)=threshP;
    } else if (pF0(i2) > (1.0-threshP)) {
      pF0(i2)=1.0-threshP;
    }
  }
  llF0=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("LLF")=LLF.head(il), Named("llF0")=llF0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}

static Eigen::VectorXd scissor_lambda_grid_cpp(int nlambda, double rlambda) {
  Eigen::VectorXd lambda(nlambda);
  if (nlambda == 1) {
    lambda(0) = 1.0;
    return lambda;
  }
  for (int i = 0; i < nlambda; ++i) {
    lambda(i) = std::pow(rlambda, static_cast<double>(i) / (nlambda - 1));
  }
  return lambda;
}

static Eigen::SparseMatrix<double> scissor_pad_omega_cpp(const Eigen::SparseMatrix<double>& omega) {
  typedef Eigen::Triplet<double> T;
  std::vector<T> triplets;
  triplets.reserve(omega.nonZeros());
  for (int col = 0; col < omega.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(omega, col); it; ++it) {
      if (it.value() != 0.0) {
        triplets.push_back(T(it.row() + 1, it.col() + 1, std::abs(it.value())));
      }
    }
  }
  Eigen::SparseMatrix<double> padded(omega.rows() + 1, omega.cols() + 1);
  padded.setFromTriplets(triplets.begin(), triplets.end());
  padded.makeCompressed();
  return padded;
}

static Eigen::MatrixXd scissor_add_intercept_cpp(const Eigen::MatrixXd& x) {
  Eigen::MatrixXd x1(x.rows(), x.cols() + 1);
  x1.col(0).setOnes();
  x1.block(0, 1, x.rows(), x.cols()) = x;
  return x1;
}

static Eigen::VectorXi scissor_col_nzero_cpp(const Eigen::MatrixXd& beta) {
  Eigen::VectorXi nzero(beta.cols());
  for (int j = 0; j < beta.cols(); ++j) {
    int n = 0;
    for (int i = 0; i < beta.rows(); ++i) {
      if (beta(i, j) != 0.0) ++n;
    }
    nzero(j) = n;
  }
  return nzero;
}

static double scissor_weighted_mean_na_cpp(const Eigen::VectorXd& x, const Eigen::VectorXd& w) {
  double num = 0.0, den = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    if (!R_IsNA(x(i)) && !R_IsNaN(x(i))) {
      num += x(i) * w(i);
      den += w(i);
    }
  }
  return den > 0.0 ? num / den : NA_REAL;
}

static int scissor_count_finite_cpp(const Eigen::VectorXd& x) {
  int n = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (!R_IsNA(x(i)) && !R_IsNaN(x(i))) ++n;
  }
  return n;
}

static Eigen::VectorXd scissor_cv_mean_cpp(const Eigen::MatrixXd& cv_rss,
                                           const Eigen::VectorXd& weights) {
  Eigen::VectorXd cvm(cv_rss.cols());
  for (int j = 0; j < cv_rss.cols(); ++j) {
    Eigen::VectorXd cvraw(cv_rss.rows());
    for (int i = 0; i < cv_rss.rows(); ++i) {
      cvraw(i) = (R_IsNA(cv_rss(i, j)) || R_IsNaN(cv_rss(i, j))) ?
        NA_REAL : cv_rss(i, j) / weights(i);
    }
    cvm(j) = scissor_weighted_mean_na_cpp(cvraw, weights);
  }
  return cvm;
}

static Eigen::VectorXd scissor_cv_se_cpp(const Eigen::MatrixXd& cv_rss,
                                         const Eigen::VectorXd& weights,
                                         const Eigen::VectorXd& cvm) {
  Eigen::VectorXd cvse(cv_rss.cols());
  for (int j = 0; j < cv_rss.cols(); ++j) {
    Eigen::VectorXd cvraw(cv_rss.rows());
    for (int i = 0; i < cv_rss.rows(); ++i) {
      cvraw(i) = (R_IsNA(cv_rss(i, j)) || R_IsNaN(cv_rss(i, j))) ?
        NA_REAL : cv_rss(i, j) / weights(i);
    }
    int nfoldi = scissor_count_finite_cpp(cvraw);
    if (nfoldi <= 1 || R_IsNA(cvm(j)) || R_IsNaN(cvm(j))) {
      cvse(j) = NA_REAL;
      continue;
    }
    Eigen::VectorXd sq(cvraw.size());
    for (int i = 0; i < cvraw.size(); ++i) {
      sq(i) = (R_IsNA(cvraw(i)) || R_IsNaN(cvraw(i))) ?
        NA_REAL : std::pow(cvraw(i) - cvm(j), 2.0);
    }
    cvse(j) = std::sqrt(scissor_weighted_mean_na_cpp(sq, weights) / (nfoldi - 1));
  }
  return cvse;
}

static int scissor_which_min_cpp(const Eigen::VectorXd& x) {
  int idx = -1;
  double best = R_PosInf;
  for (int i = 0; i < x.size(); ++i) {
    if (!R_IsNA(x(i)) && !R_IsNaN(x(i)) && x(i) < best) {
      best = x(i);
      idx = i;
    }
  }
  return idx;
}

static IntegerVector scissor_fold_rows_cpp(const IntegerVector& foldid, int fold, bool keep_equal) {
  std::vector<int> rows;
  rows.reserve(foldid.size());
  for (int i = 0; i < foldid.size(); ++i) {
    if ((foldid[i] == fold) == keep_equal) {
      rows.push_back(i);
    }
  }
  return wrap(rows);
}

static Eigen::MatrixXd scissor_subset_rows_cpp(const Eigen::MatrixXd& x, const IntegerVector& rows) {
  Eigen::MatrixXd out(rows.size(), x.cols());
  for (int i = 0; i < rows.size(); ++i) {
    out.row(i) = x.row(rows[i]);
  }
  return out;
}

static Eigen::VectorXd scissor_subset_vec_cpp(const Eigen::VectorXd& x, const IntegerVector& rows) {
  Eigen::VectorXd out(rows.size());
  for (int i = 0; i < rows.size(); ++i) {
    out(i) = x(rows[i]);
  }
  return out;
}

static Eigen::VectorXi scissor_rank_desc_abs_min_cpp(const Eigen::VectorXd& x) {
  std::vector<double> vals(x.size());
  for (int i = 0; i < x.size(); ++i) vals[i] = std::abs(x(i));
  std::vector<double> sorted = vals;
  std::sort(sorted.begin(), sorted.end(), std::greater<double>());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  Eigen::VectorXi rank(x.size());
  for (int i = 0; i < x.size(); ++i) {
    rank(i) = static_cast<int>(
      std::lower_bound(sorted.begin(), sorted.end(), vals[i], std::greater<double>()) -
        sorted.begin()
    ) + 1;
  }
  return rank;
}

static NumericVector scissor_numeric_from_eigen_cpp(const Eigen::VectorXd& x) {
  NumericVector out(x.size());
  for (int i = 0; i < x.size(); ++i) out[i] = x(i);
  return out;
}

static NumericMatrix scissor_matrix_from_eigen_cpp(const Eigen::MatrixXd& x) {
  NumericMatrix out(x.rows(), x.cols());
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) out(i, j) = x(i, j);
  }
  return out;
}

static DataFrame scissor_log_fit_df_cpp(const Eigen::VectorXd& lambda,
                                        const Eigen::VectorXd& ll,
                                        double ll0,
                                        const Eigen::VectorXd& cvm,
                                        const Eigen::VectorXd& cvse,
                                        const Eigen::VectorXi& nzero,
                                        int index_min) {
  int n = lambda.size();
  NumericVector pdev(n), cvm_r(n), cvse_r(n);
  IntegerVector nzero_r(n);
  CharacterVector index(n);
  for (int i = 0; i < n; ++i) {
    pdev[i] = (ll0 - ll(i)) / ll0;
    cvm_r[i] = cvm(i);
    cvse_r[i] = cvse(i);
    nzero_r[i] = std::max(nzero(i) - 1, 0);
    index[i] = (i == index_min) ? "*" : "";
  }
  return DataFrame::create(
    Named("lambda") = scissor_numeric_from_eigen_cpp(lambda),
    Named("pDev") = pdev,
    Named("cvm") = cvm_r,
    Named("cvse") = cvse_r,
    Named("nzero") = nzero_r,
    Named("index") = index,
    Named("stringsAsFactors") = false
  );
}

static DataFrame scissor_lm_fit_df_cpp(const Eigen::VectorXd& lambda,
                                       const Eigen::VectorXd& rsq,
                                       const Eigen::VectorXd& cvm,
                                       const Eigen::VectorXd& cvse,
                                       const Eigen::VectorXi& nzero,
                                       int index_min) {
  int n = lambda.size();
  NumericVector cvm_r(n), cvse_r(n);
  IntegerVector nzero_r(n);
  CharacterVector index(n);
  for (int i = 0; i < n; ++i) {
    cvm_r[i] = cvm(i);
    cvse_r[i] = cvse(i);
    nzero_r[i] = nzero(i);
    index[i] = (i == index_min) ? "*" : "";
  }
  return DataFrame::create(
    Named("lambda") = scissor_numeric_from_eigen_cpp(lambda),
    Named("rsq") = scissor_numeric_from_eigen_cpp(rsq),
    Named("cvm") = cvm_r,
    Named("cvse") = cvse_r,
    Named("nzero") = nzero_r,
    Named("index") = index,
    Named("stringsAsFactors") = false
  );
}

static Eigen::VectorXd scissor_trim_lm_cv_for_lambda_cpp(
    int il,
    const std::vector<Eigen::MatrixXd>& beta_folds,
    const std::vector<Eigen::MatrixXd>& beta_std_folds,
    const Eigen::VectorXi& nzero,
    const Eigen::MatrixXd& a0_folds,
    const Eigen::VectorXd& wbeta,
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& y,
    const IntegerVector& foldid,
    const Eigen::VectorXd& weights,
    int nfolds) {
  int p = x.cols();
  int numi = nzero(il);
  Eigen::MatrixXd beta_i(p, nfolds);
  Eigen::MatrixXd beta_std_i(p, nfolds);
  Eigen::VectorXi betao(nfolds);
  for (int f = 0; f < nfolds; ++f) {
    beta_i.col(f) = beta_folds[f].col(il);
    beta_std_i.col(f) = beta_std_folds[f].col(il);
    betao(f) = (beta_i.col(f).array() != 0.0).count();
  }

  int max_betao = betao.maxCoeff();
  int numi2 = std::min(max_betao, numi);
  int ntrim = std::max(numi2, 1);
  Eigen::MatrixXd cv_rss = Eigen::MatrixXd::Constant(nfolds, ntrim, NA_REAL);

  for (int f = 0; f < nfolds; ++f) {
    IntegerVector test_rows = scissor_fold_rows_cpp(foldid, f + 1, true);
    Eigen::MatrixXd x_test = scissor_subset_rows_cpp(x, test_rows);
    Eigen::VectorXd y_test = scissor_subset_vec_cpp(y, test_rows);
    int numj = std::min(betao(f), numi);

    if (numi2 <= 0 || numj == 0) {
      Eigen::VectorXd beta_empty(2);
      beta_empty.setZero();
      Eigen::VectorXi loc_empty(2);
      loc_empty.setZero();
      Eigen::VectorXd trim = cvTrimLmC(
        beta_empty, 0, numi2, loc_empty, x_test, y_test,
        test_rows.size(), a0_folds(f, il)
      );
      for (int k = 0; k < trim.size(); ++k) cv_rss(f, k) = trim(k);
      continue;
    }

    Eigen::VectorXd beta_j = beta_i.col(f);
    Eigen::VectorXd beta_std_rank = beta_std_i.col(f);
    double max_abs = beta_std_rank.cwiseAbs().maxCoeff();
    for (int k = 0; k < wbeta.size(); ++k) {
      if (wbeta(k) == 0.0) beta_std_rank(k) = max_abs + 1.0;
    }
    Eigen::VectorXi ranks = scissor_rank_desc_abs_min_cpp(beta_std_rank);
    std::vector<std::pair<int, int> > selected;
    for (int k = 0; k < ranks.size(); ++k) {
      if (ranks(k) <= numj) selected.push_back(std::make_pair(ranks(k), k));
    }
    std::sort(selected.begin(), selected.end());

    Eigen::VectorXd beta_sel(selected.size());
    Eigen::VectorXi loc(selected.size());
    for (int k = 0; k < static_cast<int>(selected.size()); ++k) {
      beta_sel(k) = beta_j(selected[k].second);
      loc(k) = selected[k].second;
    }
    Eigen::VectorXd trim = cvTrimLmC(
      beta_sel, numj, numi2, loc, x_test, y_test,
      test_rows.size(), a0_folds(f, il)
    );
    for (int k = 0; k < trim.size(); ++k) cv_rss(f, k) = trim(k);
  }

  return scissor_cv_mean_cpp(cv_rss, weights);
}

static Eigen::VectorXd scissor_trim_cv_for_lambda_cpp(
    int il,
    const std::vector<Eigen::MatrixXd>& beta_folds,
    const std::vector<Eigen::MatrixXd>& beta_std_folds,
    const Eigen::VectorXi& nzero,
    const Eigen::VectorXd& wbeta1,
    const Eigen::MatrixXd& x1,
    const Eigen::VectorXd& y,
    const IntegerVector& foldid,
    const Eigen::VectorXd& weights,
    int nfolds,
    double threshP) {
  int p1 = x1.cols();
  int numi = nzero(il);
  Eigen::MatrixXd beta_i(p1, nfolds);
  Eigen::MatrixXd beta_std_i(p1, nfolds);
  Eigen::VectorXi betao(nfolds);
  for (int f = 0; f < nfolds; ++f) {
    beta_i.col(f) = beta_folds[f].col(il);
    beta_std_i.col(f) = beta_std_folds[f].col(il);
    betao(f) = (beta_i.col(f).array() != 0.0).count();
  }
  int max_betao = betao.maxCoeff();
  int numi2 = std::max(std::min(max_betao, numi), 1);
  Eigen::MatrixXd cv_rss = Eigen::MatrixXd::Constant(nfolds, numi2, NA_REAL);
  for (int f = 0; f < nfolds; ++f) {
    IntegerVector test_rows = scissor_fold_rows_cpp(foldid, f + 1, true);
    int numj = std::min(betao(f), numi);
    Eigen::VectorXd beta_j = beta_i.col(f);
    Eigen::VectorXd beta_std_j = beta_std_i.col(f);
    Eigen::VectorXd beta_std_rank = beta_std_j;
    double max_abs = beta_std_j.cwiseAbs().maxCoeff();
    for (int k = 0; k < wbeta1.size(); ++k) {
      if (wbeta1(k) == 0.0) beta_std_rank(k) = max_abs + 1.0;
    }
    Eigen::VectorXi ranks = scissor_rank_desc_abs_min_cpp(beta_std_rank);
    std::vector<std::pair<int, int> > selected;
    for (int k = 0; k < ranks.size(); ++k) {
      if (ranks(k) <= numj) selected.push_back(std::make_pair(ranks(k), k));
    }
    std::sort(selected.begin(), selected.end());
    Eigen::VectorXd beta_sel(selected.size());
    Eigen::VectorXi loc(selected.size());
    for (int k = 0; k < static_cast<int>(selected.size()); ++k) {
      beta_sel(k) = beta_j(selected[k].second);
      loc(k) = selected[k].second;
    }
    Eigen::MatrixXd x_test = scissor_subset_rows_cpp(x1, test_rows);
    Eigen::VectorXd y_test = scissor_subset_vec_cpp(y, test_rows);
    Eigen::VectorXd trim = cvTrimLogC(beta_sel, numj, numi2, loc, x_test, y_test, test_rows.size(), threshP);
    for (int k = 0; k < trim.size(); ++k) cv_rss(f, k) = trim(k);
  }
  return scissor_cv_mean_cpp(cv_rss, weights);
}

// [[Rcpp::export]]
List scissor_gaussian_net_fit_cpp(Eigen::MatrixXd x,
                                  Eigen::VectorXd y,
                                  Eigen::SparseMatrix<double> omega,
                                  double alpha,
                                  Nullable<NumericVector> lambda = R_NilValue,
                                  int nlambda = 100,
                                  Nullable<IntegerVector> foldid = R_NilValue,
                                  bool inzero = true,
                                  bool isd = false,
                                  double thresh = 1e-7,
                                  int maxit = 100000,
                                  double threshP = 1e-5) {
  int N0 = x.rows();
  int p = x.cols();
  Eigen::MatrixXd x_original = x;
  Eigen::VectorXd y_original = y;

  bool has_lambda = lambda.isNotNull();
  Eigen::VectorXd lambda_vec;
  int ilambda = 1;
  if (has_lambda) {
    NumericVector lambda_r(lambda);
    nlambda = lambda_r.size();
    lambda_vec.resize(nlambda);
    for (int i = 0; i < nlambda; ++i) lambda_vec(i) = lambda_r[i];
    ilambda = 0;
  } else {
    double rlambda = (N0 > p) ? 0.0001 : 0.01;
    lambda_vec = scissor_lambda_grid_cpp(nlambda, rlambda);
  }

  Eigen::SparseMatrix<double> omega1 = scissor_pad_omega_cpp(omega);
  Eigen::VectorXi sgn1 = Eigen::VectorXi::Ones(p + 1);
  List W = OmegaSC(omega1, sgn1);
  Eigen::SparseMatrix<double> Womega = as<Eigen::SparseMatrix<double> >(W["Omega"]);
  Eigen::MatrixXd Wloc = as<Eigen::MatrixXd>(W["loc"]);
  Wloc.array() += 1.0;
  Eigen::VectorXi Wnadj = as<Eigen::VectorXi>(W["nadj"]);

  Eigen::VectorXd wbeta = Eigen::VectorXd::Ones(p);
  List out = NetLmC(
    x, y, alpha, lambda_vec, nlambda, ilambda, wbeta,
    Womega, Wloc, Wnadj, p, N0, thresh, maxit, threshP
  );
  int nlambdai = as<int>(out["nlambda"]);
  if (nlambdai == 0) {
    return R_NilValue;
  }

  Eigen::MatrixXd beta = as<Eigen::MatrixXd>(out["Beta"]);
  Eigen::MatrixXd beta_std = as<Eigen::MatrixXd>(out["BetaSTD"]);
  Eigen::VectorXd lambda_all = as<Eigen::VectorXd>(out["lambda"]);
  Eigen::VectorXd lambdai = lambda_all.head(nlambdai);
  Eigen::VectorXd rsq_all = as<Eigen::VectorXd>(out["rsq"]);
  Eigen::VectorXd rsq = rsq_all.head(nlambdai);
  Eigen::VectorXi flag = as<Eigen::VectorXi>(out["flag"]);
  Eigen::VectorXi nzero = scissor_col_nzero_cpp(beta.leftCols(nlambdai));

  if (foldid.isNull()) {
    Eigen::MatrixXd beta_return = isd ? beta_std.leftCols(nlambdai) : beta.leftCols(nlambdai);
    DataFrame fit = DataFrame::create(
      Named("lambda") = scissor_numeric_from_eigen_cpp(lambdai),
      Named("rsq") = scissor_numeric_from_eigen_cpp(rsq),
      Named("nzero") = nzero,
      Named("stringsAsFactors") = false
    );
    return List::create(
      Named("Beta") = scissor_matrix_from_eigen_cpp(beta_return),
      Named("fit") = fit,
      Named("penalty") = "Net",
      Named("adaptive") = LogicalVector::create(false, false),
      Named("flag") = flag
    );
  }

  IntegerVector foldid_r(foldid);
  int nfolds = 0;
  for (int i = 0; i < foldid_r.size(); ++i) nfolds = std::max(nfolds, foldid_r[i]);
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(nfolds);
  for (int i = 0; i < foldid_r.size(); ++i) weights(foldid_r[i] - 1) += 1.0;

  std::vector<Eigen::MatrixXd> beta_folds(nfolds, Eigen::MatrixXd::Zero(p, nlambdai));
  std::vector<Eigen::MatrixXd> beta_std_folds(nfolds, Eigen::MatrixXd::Zero(p, nlambdai));
  Eigen::MatrixXd a0_folds = Eigen::MatrixXd::Zero(nfolds, nlambdai);
  Eigen::MatrixXd cv_rss = Eigen::MatrixXd::Constant(nfolds, nlambdai, NA_REAL);

  for (int f = 0; f < nfolds; ++f) {
    IntegerVector train_rows = scissor_fold_rows_cpp(foldid_r, f + 1, false);
    IntegerVector test_rows = scissor_fold_rows_cpp(foldid_r, f + 1, true);
    Eigen::MatrixXd x_train = scissor_subset_rows_cpp(x_original, train_rows);
    Eigen::VectorXd y_train = scissor_subset_vec_cpp(y_original, train_rows);
    Eigen::MatrixXd x_test = scissor_subset_rows_cpp(x_original, test_rows);
    Eigen::VectorXd y_test = scissor_subset_vec_cpp(y_original, test_rows);

    List fold_out = cvNetLmC(
      x_train, y_train, alpha, lambdai, nlambdai, wbeta,
      Womega, Wloc, Wnadj, train_rows.size(), p, thresh, maxit,
      x_test, y_test, test_rows.size(), threshP
    );
    int fold_nlambda = as<int>(fold_out["nlambda"]);
    Eigen::MatrixXd fold_beta = as<Eigen::MatrixXd>(fold_out["Beta"]);
    Eigen::MatrixXd fold_beta_std = as<Eigen::MatrixXd>(fold_out["BetaSTD"]);
    Eigen::VectorXd fold_rssp = as<Eigen::VectorXd>(fold_out["RSSp"]);
    Eigen::VectorXd fold_a0 = as<Eigen::VectorXd>(fold_out["a0S"]);
    for (int j = 0; j < fold_nlambda; ++j) {
      beta_folds[f].col(j) = fold_beta.col(j);
      beta_std_folds[f].col(j) = fold_beta_std.col(j);
      cv_rss(f, j) = fold_rssp(j);
      a0_folds(f, j) = fold_a0(j);
    }
  }

  Eigen::VectorXd cvm = scissor_cv_mean_cpp(cv_rss, weights);
  Eigen::VectorXd cvse = scissor_cv_se_cpp(cv_rss, weights, cvm);
  int index_min = scissor_which_min_cpp(cvm);
  DataFrame fit = scissor_lm_fit_df_cpp(lambdai, rsq, cvm, cvse, nzero, index_min);

  if (!inzero) {
    Eigen::VectorXd beta_return = isd ? beta_std.col(index_min) : beta.col(index_min);
    return List::create(
      Named("Beta") = scissor_numeric_from_eigen_cpp(beta_return),
      Named("fit") = fit,
      Named("lambda.min") = lambdai(index_min),
      Named("penalty") = "Net",
      Named("adaptive") = LogicalVector::create(false, false),
      Named("flag") = flag
    );
  }

  int il0 = index_min;
  Eigen::VectorXd cv_min = Eigen::VectorXd::Constant(nlambdai, NA_REAL);
  std::vector<Eigen::VectorXd> trim_cvm(nlambdai);
  while (true) {
    trim_cvm[il0] = scissor_trim_lm_cv_for_lambda_cpp(
      il0, beta_folds, beta_std_folds, nzero, a0_folds, wbeta,
      x_original, y_original, foldid_r, weights, nfolds
    );
    Eigen::VectorXd temi = trim_cvm[il0];
    cv_min(il0) = temi.minCoeff();

    int neigh[2] = {il0 - 1, il0 + 1};
    for (int ni = 0; ni < 2; ++ni) {
      int il1 = neigh[ni];
      if (il1 < 0 || il1 >= nlambdai) break;
      if (R_IsNA(cv_min(il1)) || R_IsNaN(cv_min(il1))) {
        trim_cvm[il1] = scissor_trim_lm_cv_for_lambda_cpp(
          il1, beta_folds, beta_std_folds, nzero, a0_folds, wbeta,
          x_original, y_original, foldid_r, weights, nfolds
        );
        cv_min(il1) = trim_cvm[il1].minCoeff();
      }
    }

    int new_il0 = scissor_which_min_cpp(cv_min);
    if (new_il0 == il0) break;
    il0 = new_il0;
  }

  int index0 = scissor_which_min_cpp(cv_min);
  Eigen::VectorXd beta0 = beta.col(index0);
  Eigen::VectorXd beta_std0 = beta_std.col(index0);
  Eigen::VectorXd beta0_rank = beta_std.col(index0);
  double max_abs = beta0_rank.cwiseAbs().maxCoeff();
  for (int i = 0; i < wbeta.size(); ++i) {
    if (wbeta(i) == 0.0) beta0_rank(i) = max_abs + 1.0;
  }

  Eigen::VectorXd temi = trim_cvm[index0];
  int cuti_zero_based = scissor_which_min_cpp(temi);
  std::vector<double> abs_vals(beta0_rank.size());
  for (int i = 0; i < beta0_rank.size(); ++i) abs_vals[i] = std::abs(beta0_rank(i));
  std::sort(abs_vals.begin(), abs_vals.end(), std::greater<double>());
  int threshold_idx = std::min(cuti_zero_based + 1, static_cast<int>(abs_vals.size()) - 1);
  double threshold = abs_vals[threshold_idx];
  for (int i = 0; i < beta0.size(); ++i) {
    if (std::abs(beta0_rank(i)) <= threshold) {
      beta0(i) = 0.0;
      beta_std0(i) = 0.0;
    }
  }

  DataFrame fit0 = DataFrame::create(
    Named("lambda") = NumericVector::create(lambdai(index0)),
    Named("cvm") = NumericVector::create(cv_min(index0)),
    Named("nzero") = IntegerVector::create((beta0.array() != 0.0).count()),
    Named("stringsAsFactors") = false
  );
  Eigen::VectorXd beta_return = isd ? beta_std.col(index_min) : beta.col(index_min);
  Eigen::VectorXd beta0_return = isd ? beta_std0 : beta0;
  return List::create(
    Named("Beta") = scissor_numeric_from_eigen_cpp(beta_return),
    Named("Beta0") = scissor_numeric_from_eigen_cpp(beta0_return),
    Named("fit") = fit,
    Named("fit0") = fit0,
    Named("lambda.min") = lambdai(index_min),
    Named("lambda.opt") = lambdai(index0),
    Named("penalty") = "Net",
    Named("adaptive") = LogicalVector::create(false, false),
    Named("flag") = flag
  );
}

// [[Rcpp::export]]
List scissor_binomial_net_fit_cpp(Eigen::MatrixXd x,
                                  Eigen::VectorXd y,
                                  Eigen::SparseMatrix<double> omega,
                                  double alpha,
                                  Nullable<NumericVector> lambda = R_NilValue,
                                  int nlambda = 100,
                                  Nullable<IntegerVector> foldid = R_NilValue,
                                  bool inzero = true,
                                  bool isd = false,
                                  double thresh = 1e-7,
                                  int maxit = 100000,
                                  double threshP = 1e-5) {
  int N0 = x.rows();
  int p = x.cols();
  bool has_lambda = lambda.isNotNull();
  Eigen::VectorXd lambda_vec;
  int ilambda = 1;
  if (has_lambda) {
    NumericVector lambda_r(lambda);
    nlambda = lambda_r.size();
    lambda_vec.resize(nlambda);
    for (int i = 0; i < nlambda; ++i) lambda_vec(i) = lambda_r[i];
    ilambda = 0;
  } else {
    double rlambda = (N0 > p) ? 0.0001 : 0.01;
    lambda_vec = scissor_lambda_grid_cpp(nlambda, rlambda);
  }

  Eigen::SparseMatrix<double> omega1 = scissor_pad_omega_cpp(omega);
  Eigen::VectorXi sgn1 = Eigen::VectorXi::Ones(p + 1);
  List W = OmegaSC(omega1, sgn1);
  Eigen::SparseMatrix<double> Womega = as<Eigen::SparseMatrix<double> >(W["Omega"]);
  Eigen::MatrixXd Wloc = as<Eigen::MatrixXd>(W["loc"]);
  Wloc.array() += 1.0;
  Eigen::VectorXi Wnadj = as<Eigen::VectorXi>(W["nadj"]);

  Eigen::MatrixXd x1 = scissor_add_intercept_cpp(x);
  int p1 = p + 1;
  Eigen::ArrayXd wbeta1 = Eigen::ArrayXd::Ones(p1);
  Eigen::ArrayXd wbetai = Eigen::ArrayXd::Ones(p1);
  wbeta1(0) = 0.0;
  wbetai(0) = 0.0;
  Eigen::VectorXd wbeta1_vec = wbeta1.matrix();

  List out = NetLogC(
    x1, y, alpha, lambda_vec, nlambda, ilambda, wbeta1, wbetai,
    Womega, Wloc, Wnadj, p1, N0, thresh, maxit, threshP
  );
  int nlambdai = as<int>(out["nlambda"]);
  if (nlambdai == 0) {
    return R_NilValue;
  }
  Eigen::MatrixXd beta = as<Eigen::MatrixXd>(out["Beta"]);
  Eigen::MatrixXd beta_std = as<Eigen::MatrixXd>(out["BetaSTD"]);
  Eigen::VectorXd lambda_all = as<Eigen::VectorXd>(out["lambda"]);
  Eigen::VectorXd lambdai = lambda_all.head(nlambdai);
  Eigen::VectorXd ll = as<Eigen::VectorXd>(out["LL"]);
  double ll0 = as<double>(out["ll0"]);
  Eigen::VectorXi flag = as<Eigen::VectorXi>(out["flag"]);
  Eigen::VectorXi nzero = scissor_col_nzero_cpp(beta.leftCols(nlambdai));

  if (foldid.isNull()) {
    Eigen::MatrixXd beta_return = isd ? beta_std.leftCols(nlambdai) : beta.leftCols(nlambdai);
    NumericVector pdev(nlambdai);
    IntegerVector nzero_r(nlambdai);
    for (int i = 0; i < nlambdai; ++i) {
      pdev[i] = (ll0 - ll(i)) / ll0;
      nzero_r[i] = std::max(nzero(i) - 1, 0);
    }
    DataFrame fit = DataFrame::create(
      Named("lambda") = scissor_numeric_from_eigen_cpp(lambdai),
      Named("pDev") = pdev,
      Named("nzero") = nzero_r,
      Named("stringsAsFactors") = false
    );
    return List::create(
      Named("Beta") = scissor_matrix_from_eigen_cpp(beta_return),
      Named("fit") = fit,
      Named("penalty") = "Net",
      Named("adaptive") = LogicalVector::create(false, false),
      Named("flag") = flag
    );
  }

  IntegerVector foldid_r(foldid);
  int nfolds = 0;
  for (int i = 0; i < foldid_r.size(); ++i) nfolds = std::max(nfolds, foldid_r[i]);
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(nfolds);
  for (int i = 0; i < foldid_r.size(); ++i) weights(foldid_r[i] - 1) += 1.0;

  std::vector<Eigen::MatrixXd> beta_folds(nfolds, Eigen::MatrixXd::Zero(p1, nlambdai));
  std::vector<Eigen::MatrixXd> beta_std_folds(nfolds, Eigen::MatrixXd::Zero(p1, nlambdai));
  Eigen::MatrixXd cv_rss = Eigen::MatrixXd::Constant(nfolds, nlambdai, NA_REAL);

  for (int f = 0; f < nfolds; ++f) {
    IntegerVector train_rows = scissor_fold_rows_cpp(foldid_r, f + 1, false);
    IntegerVector test_rows = scissor_fold_rows_cpp(foldid_r, f + 1, true);
    Eigen::VectorXd y_train = scissor_subset_vec_cpp(y, train_rows);
    bool has0 = false, has1 = false;
    for (int i = 0; i < y_train.size(); ++i) {
      has0 = has0 || y_train(i) == 0.0;
      has1 = has1 || y_train(i) == 1.0;
    }
    if (!has0 || !has1) continue;

    Eigen::MatrixXd x_train = scissor_subset_rows_cpp(x1, train_rows);
    Eigen::MatrixXd x_test = scissor_subset_rows_cpp(x1, test_rows);
    Eigen::VectorXd y_test = scissor_subset_vec_cpp(y, test_rows);
    Eigen::ArrayXd wbeta_fold = wbeta1;
    Eigen::ArrayXd wbetai_fold = wbetai;
    List fold_out = cvNetLogC(
      x_train, y_train, alpha, lambdai, nlambdai, wbeta_fold, wbetai_fold,
      Womega, Wloc, Wnadj, p1, train_rows.size(), thresh, maxit,
      x_test, y_test, test_rows.size(), threshP
    );
    int fold_nlambda = as<int>(fold_out["nlambda"]);
    Eigen::MatrixXd fold_beta = as<Eigen::MatrixXd>(fold_out["Beta"]);
    Eigen::MatrixXd fold_beta_std = as<Eigen::MatrixXd>(fold_out["BetaSTD"]);
    Eigen::VectorXd fold_llf = as<Eigen::VectorXd>(fold_out["LLF"]);
    for (int j = 0; j < fold_nlambda; ++j) {
      beta_folds[f].col(j) = fold_beta.col(j);
      beta_std_folds[f].col(j) = fold_beta_std.col(j);
      cv_rss(f, j) = 2.0 * (0.0 - fold_llf(j)) * test_rows.size();
    }
  }

  Eigen::VectorXd cvm = scissor_cv_mean_cpp(cv_rss, weights);
  Eigen::VectorXd cvse = scissor_cv_se_cpp(cv_rss, weights, cvm);
  int index_min = scissor_which_min_cpp(cvm);
  DataFrame fit = scissor_log_fit_df_cpp(
    lambdai, ll, ll0, cvm, cvse, nzero, index_min
  );

  if (!inzero) {
    Eigen::MatrixXd beta_return = isd ? beta_std.col(index_min) : beta.col(index_min);
    return List::create(
      Named("Beta") = scissor_matrix_from_eigen_cpp(beta_return),
      Named("fit") = fit,
      Named("lambda.min") = lambdai(index_min),
      Named("penalty") = "Net",
      Named("adaptive") = LogicalVector::create(false, false),
      Named("flag") = flag
    );
  }

  int il0 = index_min;
  Eigen::VectorXd cv_min = Eigen::VectorXd::Constant(nlambdai, NA_REAL);
  std::vector<Eigen::VectorXd> trim_cvm(nlambdai);
  while (true) {
    trim_cvm[il0] = scissor_trim_cv_for_lambda_cpp(
      il0, beta_folds, beta_std_folds, nzero, wbeta1_vec, x1, y,
      foldid_r, weights, nfolds, threshP
    );
    Eigen::VectorXd temi = trim_cvm[il0];
    if (temi.size() > 1) {
      double best = R_PosInf;
      for (int i = 1; i < temi.size(); ++i) {
        if (!R_IsNA(temi(i)) && !R_IsNaN(temi(i)) && temi(i) < best) best = temi(i);
      }
      cv_min(il0) = best;
    } else {
      cv_min(il0) = temi(0);
    }

    int neigh[2] = {il0 - 1, il0 + 1};
    for (int ni = 0; ni < 2; ++ni) {
      int il1 = neigh[ni];
      if (il1 < 0 || il1 >= nlambdai) break;
      if (R_IsNA(cv_min(il1)) || R_IsNaN(cv_min(il1))) {
        trim_cvm[il1] = scissor_trim_cv_for_lambda_cpp(
          il1, beta_folds, beta_std_folds, nzero, wbeta1_vec, x1, y,
          foldid_r, weights, nfolds, threshP
        );
        Eigen::VectorXd temi1 = trim_cvm[il1];
        if (temi1.size() > 1) {
          double best = R_PosInf;
          for (int i = 1; i < temi1.size(); ++i) {
            if (!R_IsNA(temi1(i)) && !R_IsNaN(temi1(i)) && temi1(i) < best) best = temi1(i);
          }
          cv_min(il1) = best;
        } else {
          cv_min(il1) = temi1(0);
        }
      }
    }
    int new_il0 = scissor_which_min_cpp(cv_min);
    if (new_il0 == il0) break;
    il0 = new_il0;
  }

  int index0 = scissor_which_min_cpp(cv_min);
  Eigen::VectorXd beta0 = beta.col(index0);
  Eigen::VectorXd beta_std0 = beta_std.col(index0);
  Eigen::VectorXd temi = trim_cvm[index0];
  int cuti_one_based = 1;
  if (temi.size() > 1) {
    double best = R_PosInf;
    int best_idx = 1;
    for (int i = 1; i < temi.size(); ++i) {
      if (!R_IsNA(temi(i)) && !R_IsNaN(temi(i)) && temi(i) < best) {
        best = temi(i);
        best_idx = i;
      }
    }
    cuti_one_based = best_idx + 1;
  }

  Eigen::VectorXd beta0_rank = beta_std.col(index0);
  double max_abs = beta0_rank.cwiseAbs().maxCoeff();
  for (int i = 0; i < wbeta1_vec.size(); ++i) {
    if (wbeta1_vec(i) == 0.0) beta0_rank(i) = max_abs + 1.0;
  }
  std::vector<double> abs_vals(beta0_rank.size());
  for (int i = 0; i < beta0_rank.size(); ++i) abs_vals[i] = std::abs(beta0_rank(i));
  std::sort(abs_vals.begin(), abs_vals.end(), std::greater<double>());
  int threshold_idx = std::min(cuti_one_based, static_cast<int>(abs_vals.size()) - 1);
  double threshold = abs_vals[threshold_idx];
  for (int i = 0; i < beta0.size(); ++i) {
    if (std::abs(beta0_rank(i)) <= threshold) {
      beta0(i) = 0.0;
      beta_std0(i) = 0.0;
    }
  }

  DataFrame fit0 = DataFrame::create(
    Named("lambda") = NumericVector::create(lambdai(index0)),
    Named("cvm") = NumericVector::create(cv_min(index0)),
    Named("nzero") = IntegerVector::create(cuti_one_based - 1),
    Named("stringsAsFactors") = false
  );
  Eigen::MatrixXd beta_return = isd ? beta_std.col(index_min) : beta.col(index_min);
  Eigen::VectorXd beta0_return = isd ? beta_std0 : beta0;
  return List::create(
    Named("Beta") = scissor_matrix_from_eigen_cpp(beta_return),
    Named("Beta0") = scissor_numeric_from_eigen_cpp(beta0_return),
    Named("fit") = fit,
    Named("fit0") = fit0,
    Named("lambda.min") = lambdai(index_min),
    Named("lambda.opt") = lambdai(index0),
    Named("penalty") = "Net",
    Named("adaptive") = LogicalVector::create(false, false),
    Named("flag") = flag
  );
}
