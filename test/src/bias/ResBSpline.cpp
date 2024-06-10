/*
    created: July 21, 2017 
    last modified: July 21, 2017
    author: Chenghan Li (lch004218@gmail.com)

    Enable plumed to use BSplines to compute bias potential and thus bias force
    U_\text{bias}(cv)=scaling_factor*\sum_{i=1} coeffs[i]*basis_function[i](cv)
    where scaling_factor=coeffs[0]
*/
#include "Bias.h"
#include "ActionRegister.h"
#include "tools/BSpline.h"

using namespace std;

namespace PLMD {
namespace bias {

class ResBSpline : public Bias {
    int order,num_bins;//default: order=6 num_bins=25
    double cut_lo,cut_hi;
    std::vector<double> coeffs;
    BSpline* table;
    Value* valueForce2;
public:
    explicit ResBSpline(const ActionOptions&);
    ~ResBSpline();
    void calculate();
    static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ResBSpline,"RESBSPLINE")

void ResBSpline::registerKeywords(Keywords& keys) {
    Bias::registerKeywords(keys);
    keys.use("ARG");
    keys.add("compulsory","CUTLO","0.0","sets the lower cutoff of the splines");
    keys.add("compulsory","CUTHI","0.0","sets the higher cutoff of the splines");
    keys.add("optional","ORDER","specifies the order of splines");
    keys.add("optional","NUMBINS","specifies the number of bins used");
    keys.add("compulsory","COEFFS","specifies the coefficients of ths splines; the first number is the scaling factor");
    keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}

ResBSpline::ResBSpline(const ActionOptions& ao):
    PLUMED_BIAS_INIT(ao),
    cut_lo(0.0),
    cut_hi(0.0),
    order(6),
    num_bins(25)
{
    //parse
    parse("CUTLO",cut_lo);
    parse("CUTHI",cut_hi);
    parse("ORDER",order);
    parse("NUMBINS",num_bins);
    parseVector("COEFFS",coeffs);
    checkRead();

    //error
    if(cut_lo>=cut_hi) error("The outer cutoff should be larger than the inner cutoff");
    if(order<=0) error("The order of splines should be larger than 0");
    if(num_bins<=0) error("The number of bins should be larger than 0");
    if(coeffs.size()==0) error("The coefficients of splines need to be specified\n");
    if(coeffs.size()!=order+num_bins+1) error("Inconsistent number of coefficients with number of bins and order");


    table=new BSpline(order,cut_lo,cut_hi,num_bins);
    //write log
    log.printf("  using %d basis functions of order %d on %d bins from %lf to %lf\n",order+num_bins,order,num_bins,cut_lo,cut_hi);
    log.printf("  scaling factor=\t%lf\n",coeffs[0]);
    for(int i=1;i<=order+num_bins;++i) log.printf("  coeff%d=\t%lf\n",i,coeffs[i]);

    addComponent("force2");
    componentIsNotPeriodic("force2");
    valueForce2=getPntrToComponent("force2");
}

ResBSpline::~ResBSpline() {
    delete table;
}

void ResBSpline::calculate() {
    double ene=0.0,fcv=0.0;

//printf("v0\n");

    //table->norm=0.0;
    double bsplinenorm=table->Compute_Basis(getArgument(0));
//    if(fabs(getArgument(0)-3.913736)<1e-5) {
//        for(int j=1;j<order+num_bins;++j) printf("%lf\n",table->basis[j-1]);
//        printf("\n");
//    }
   //     for(int j=1;j<=order*2+num_bins+1;++j) { 
   //         log.printf("%d\t%lf\n",j-1,table->basis[j-1]);
   //     }
   //     for(int j=1;j<=order*2+num_bins+1;++j) { 
   //         log.printf("\t%d\t%lf\n",j-1,table->basis_deriv[j-1]);
   //     }
    if(fabs(bsplinenorm-1.0)>1e-8) { 
        char tmp[100];
        sprintf(tmp,"B-Spline normalization error cv=%lf cutoffs=%lf,%lf norm=%lf",getArgument(0),cut_lo,cut_hi,bsplinenorm);
        for(int j=1;j<=order*2+num_bins+1;++j) { 
            printf("%lf\n",table->basis[j-1]);
        }
        plumed_merror(tmp);
    }
    for(int j=1;j<=order+num_bins;++j) { 
        ene+=coeffs[j]*table->basis[j-1];
        fcv+=coeffs[j]*table->basis_deriv[j-1];
    }
    ene*= coeffs[0]; fcv*=-coeffs[0];
    setBias(ene); setOutputForce(0,fcv); valueForce2->set(fcv*fcv);
    //for debug
    //printf("raw ene=%lf kj/mol, raw fcv=%lf kj/mol/nm\n",ene,fcv);
}

}
}
