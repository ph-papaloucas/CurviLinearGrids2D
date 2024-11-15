#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
/**
    READ ME-----READ ME-----READ ME-----READ ME-----READ ME-----READ ME

    For the script to run correctly, you need to have in the same path a folder named "output"
    Then you can find the initial and final grid in that folder.


    The user has to specify in the main{} these things:
    1.//Cosine tube characteristics
        everything is self explanatory in the script
    2.//Characteristics of tuning
        ..tuner = grid.TunerSorenson1 for Thomas & Middlecoff, TunerSorenson2 for Sorenson
        ..tunerf1 and tunerf2 are the coefficients that affect the inner nodes. IF you choose 1,
                then its close to a Laplacian grid
        the following affect only TunerSOrenson2
        ..frelax_s  the initial relaxation in calculating xetaeta yetaeta (needs to be small e.g. 0.001)
        ..frelax_f the final relaxation of the above (eg 0.4)
        ..sorenson_multiplier: it affects the rate at which frelax converges to frelax_f (small = slower converges)
        if method fails to converges, lower everything of the above
    3.//Solver Characteristics
        ..solver = grid.SolverImplicit; for ADI solver, or grid.PointJacobi for point jacobi
        the rest are self explanatory

**/

///
//+++++++++++++++++++++++++++++++++++++++++++++++++ OTHER FUNCTIONS     +++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++ OTHER FUNCTIONS     +++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++ OTHER FUNCTIONS     +++++++++++++++++++++++++++++++++++++++++++++++++
///
double trap(std::vector<double>& x, std::vector<double>& y) {
    double sum=0.e0;
	for(int i=1;i<x.size();i++)
	{
        	if(x[i]<x[i-1])
		{
		  std::cout<<"x NOT in INCREASING order"<< std::endl;
		  break;
		}
        	sum=sum+.5e0*(x[i]-x[i-1])*(y[i]+y[i-1]);
      	}
      	return sum;
}

void tridiagonal_solver(const std::vector<double>& aa, //subdiagonal (lower) size n-1    **OR size n but put 0 on first element**
                        const std::vector<double>& b, //main diagonal size n
                        const std::vector<double>& cc, //superdiagonal (upper) size n-1  **OR size n but put 0 on last element**
                        const std::vector<double>& d,// RHS size n
                        std::vector<double>& ef){ //solution we are trying to compute, size n
    size_t n = ef.size(); // Get the size of the vectors b,d,f
    std::vector<double> a,c;
    a=aa;
    c=cc;

    if (aa.size() == n - 1) {     /// check aa and cc matrices if they are n-1 (as they should), or n(because its easier for the input of the user)
        a = aa;
    } else if (aa.size() == n && a[0]==0  ){
        // If aa has n elements, copy all but the first element to 'a'
        a.assign(aa.begin() + 1, aa.end());
    } else {
        std::cout <<" Tridiag ERROR: lower diagonal vector isnt given correctly " << std::endl;
        std::cout <<" You should either give it as n-1, or n with a[0]=0" << std::endl;
    }

    if (cc.size() == n - 1) {    // Ensure that 'c' has n-1 elements
        c = cc;
    } else if ( cc.size() == n && cc[n-1] == 0) {        // If cc has n elements, copy all but the last element to 'c'
        c.assign(cc.begin(), cc.end() - 1);
    } else {
        std::cout <<" Tridiag ERROR: upper diagonal vector isnt given correctly " << std::endl;
        std::cout <<" You should either give it as n-1, or n with c[n-1]=0" << std::endl;
    }
    /// check if a b c have the correct dimensions
    std::vector<double> c_star(n-1,0);    //create temporary c, d vectors so i dont have to change the others
    std::vector<double> d_star(n,0);    // it also initializes them with zeros

    c_star[0] = c[0] / b[0]; //step 2, initialize c and d vectors
    d_star[0] = d[0] / b[0];

    for (int i=1; i<n;i++){ //step2, forward sweep (Forward elimination)
        double m=b[i]-a[i-1]*c_star[i-1];
        if (i<n-1){
            c_star[i]=c[i]/m;
        }
        d_star[i]=(d[i]-a[i-1]*d_star[i-1])/m;
    }

    ef[n-1]= d_star[n-1];//step3, reverse sweep (Backwards substitution), used to update the solution vector f
    for (int i=n-2;i>=0;i--){
        ef[i]=d_star[i] - c_star[i]*ef[i+1];
    }

}//end of tridiagonal_solver

///
//+++++++++++++++++++++++++++++++++++++++++++++++++ MAIN CLASS  : GRid2DCosineTube   +++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++ MAIN CLASS  : GRid2DCosineTube   +++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++ MAIN CLASS  : GRid2DCosineTube   +++++++++++++++++++++++++++++++++++++++++++++++++
///
class GRid2DCosineTube{
public:
    std::vector<std::vector<double>> x, y;          //position
    std::vector<std::vector<double>> A, B, C, J;    //metrics g22, g12, g11,  J
    std::vector<std::vector<double>> f1, f2;        //weights (they tune the grid)
    std::vector<double> xetaeta, yetaeta;// we want these because we want relaxation in Sorenson2 ( so we need to know the previous value
    int IM,JM;                                      //size
    double residual, surface, f_residual, f_residual2;

    const int SolverPointJacobi = 0;
    const int SolverImplicit = 1;
    const int TunerSorenson1 = 1;
    const int TunerSorenson2 = 2;

    int kavalaei = 0;

private:
    int p_tuner = 0;
    double p_s_percent, p_tuner_f1, p_tuner_f2, p_s;
    double sorenson_relax; double sorenson_relax_final; double sorenson_multiplier; int sorenson_counter = 0;

    double a0, a1,d_tube, l_tube, l_neck, d_neck; // for the cosine tube(see CreateCosineTube function)
    int idx_start; //index at which the neck starts to form


public:
    GRid2DCosineTube(int IM, int JM, double L, double D, double l_of_neck, double d_of_neck) : ///Constructor
        x(IM, std::vector<double>(JM))  ,   y(IM, std::vector<double>(JM)),
        A(IM, std::vector<double>(JM))  ,   B(IM, std::vector<double>(JM)),
        C(IM, std::vector<double>(JM))  ,   J(IM, std::vector<double>(JM)),
        f1(IM, std::vector<double>(JM)) ,   f2(IM, std::vector<double>(JM)),
        xetaeta(IM), yetaeta(IM),
        IM(IM), JM(JM), residual(100.0) ,   f_residual(100.0), sorenson_relax(0.000001),
        a0(0), a1(0), l_neck(l_of_neck) ,   d_neck(d_of_neck), l_tube(L), d_tube(D), surface(0)
        {
            CreateCosineTube();
        }

    void Print2FileXY(const std::string& filename) {
        std::string outputDirectory = "output/";
        std::string filename_x = outputDirectory + filename + "_x.txt";
        std::string filename_y = outputDirectory + filename + "_y.txt";
        std::ofstream outputFileX(filename_x);
        std::ofstream outputFileY(filename_y);
        if (!outputFileX.is_open() || !outputFileY.is_open()) {
            std::cerr << "Error opening files: " << filename_x << " or " << filename_y << std::endl;
            return;
        }

        // Set precision for output
        outputFileX << std::fixed << std::setprecision(6);
        outputFileY << std::fixed << std::setprecision(6);

        // Iterate through the grid and write x and y values to the files
        for (int j = JM-1; j >= 0 ; j--) {
            for (int i = 0; i < IM; i++) {
                outputFileX << x[i][j] << " ";
                outputFileY << y[i][j] << " ";
            }
            outputFileX << std::endl;
            outputFileY << std::endl;
        }
        std::cout << "Grid has been printed to " << filename << "_ .txt " << std::endl;
    }//end of print2file
    void Print2FileF1F2(const std::string& filename) {
        std::string outputDirectory = "output/";
        std::string filename_f1 = outputDirectory + filename + "_f1.txt";
        std::string filename_f2 = outputDirectory + filename + "_f2.txt";
        std::ofstream outputFileF1(filename_f1);
        std::ofstream outputFileF2(filename_f2);
        if (!outputFileF1.is_open() || !outputFileF2.is_open()) {
            std::cerr << "Error opening files: " << filename_f1 << " or " << filename_f2 << std::endl;
            return;
        }
        // Set precision for output
        outputFileF1 << std::fixed << std::setprecision(6);
        outputFileF2 << std::fixed << std::setprecision(6);
        // Iterate through the grid and write x and y values to the files
        for (int j = JM-1; j >= 0 ; j--) {
            for (int i = 0; i < IM; i++) {
                outputFileF1 << f1[i][j] << " ";
                outputFileF2 << f2[i][j] << " ";
            }
            outputFileF1 << std::endl;
            outputFileF2 << std::endl;
        }
        std::cout << "Grid's f1 & f1 have been printed to " << filename << "_ .txt " << std::endl;
    }//end of print2file
    void Print2FileJ(const std::string& filename) {
            std::string outputDirectory = "output/";
            std::string filename_J = outputDirectory + filename + "_J.txt";
            std::ofstream outputFileJ(filename_J);
            if (!outputFileJ.is_open() ) {
                std::cerr << "Error opening files: " << filename_J << std::endl;
                return;
            }
            // Set precision for output
            outputFileJ << std::fixed << std::setprecision(6);
            // Iterate through the grid and write x and y values to the files

            for (int j = JM-1; j >= 0 ; j--) {
                for (int i = 0; i < IM; i++) {
                    outputFileJ << J[i][j] << " ";
                }
                outputFileJ << std::endl;
            }
            std::cout << "Grid's J has been printed to " << filename << "_ .txt " << std::endl;
        }//end of print2file

    void ComputeF(int TuningMethod, const double& s_percent, const double& tuner_f1, const double& tuner_f2,
                  double frelax_s, double frelax_f, double sorenson_multiplier_coeff){
        sorenson_multiplier = sorenson_multiplier_coeff;
        p_tuner = TuningMethod;
        p_s_percent = s_percent; p_tuner_f1 = tuner_f1; p_tuner_f2 = tuner_f2;
        sorenson_relax = frelax_s; sorenson_relax_final= frelax_f;
        SorensonComputeF(TunerSorenson1); //initialize
        if (TuningMethod == 1){
            sorenson_relax = sorenson_relax_final;
            std::cout << "Method: Thomas & Middlecoff" << std::endl;
        }
        if ( (p_tuner != TunerSorenson1) && (p_tuner != TunerSorenson2) ){
            std::cout << "ERROR IN GRid2D.ComputeF: Method Doesnt exists " << std::endl;
        }
        Print2FileF1F2("init");
    }

    void solve(int solver, const double res_tol, const double max_iterations,
                            const double relaxation, const int info_counter ){

        if (solver == SolverPointJacobi){
            PoisonGridPointJacobi(res_tol, max_iterations, relaxation, info_counter);
        }else if (solver == SolverImplicit){

            PoisonGridImplicit(res_tol, max_iterations, relaxation, info_counter);
        }else{
            std::cout << "ERROR IN GRid2D.solve: Solver does't exist " << std::endl;
        }
    }

    void PostProcess(){
        SurfaceJacobean();
        Print2FileF1F2("final");
        double surface2 = SurfaceTrap();
        std::cout << std::endl;
        std::cout << "Post Process REPORT " << std::endl;
        if (kavalaei == 1){
            std::cerr << " ERRORRRR: TO PLEGMA KAVALAEIIIII" << std::endl;
        }
        std::cout << "--> Jacobean surface: " << surface << std::endl;
        std::cout << "--> Difference in surfaces with trap method: " << abs(surface-surface2) << std::endl;
        std::cout << "--> Solver residual: " << residual << std::endl;
        //std::cout << "--> f_residual (change from prev. iteration) " << f_residual << std::endl;
        //std::cout << "    (if sorenson 1 was chosen, then f_residual will be 100 by default)" << std::endl;
        std::cout << "--> Sorenson relaxation : " << sorenson_relax << std::endl;
        std::cout << "--> Sorenson final relaxation : " << sorenson_relax_final << std::endl;
        std::cout << "--> Iterations to reach final relaxation : " << sorenson_counter << std::endl;

        //print2file f1,f2, s

        std::ofstream outputFile("output/tuner.txt");

        // Check if the file is open
        if (outputFile.is_open()) {
            // Write each double to a separate row
            outputFile << p_tuner_f1 << std::endl;
            outputFile << p_tuner_f2 << std::endl;
            outputFile << p_s << std::endl;

            // Close the file
            outputFile.close();

            std::cout << "Data has been written to the file." << std::endl;
        } else {
            std::cerr << "Unable to open the file." << std::endl;
        }




    }//end of PostProcess

private:
    void CreateCosineTube(){
        a0 = (d_tube-d_neck)/4;     // for cosine (y = a0*cos(a1*...))
        a1 = 2*M_PI/l_neck;    // for cosine
        const double d = (l_tube-l_neck)/2;      // so L = l_neck + d + d
        const double dx = l_tube/(IM-1); const double dy = d_tube/(JM-1);
        double x_value,y_value;
        for (int i=0; i<IM;i++){ //create boundary nodes top/bot
            x_value = dx*i;
            x[i][0]    = x_value;
            x[i][JM-1] = x_value;
            if ( (x_value > d) && (x_value < (d + l_neck))  ){
                y[i][0]     =   - a0*cos(a1*(x_value - d)) + a0;
                y[i][JM-1]  =   + a0*cos(a1*(x_value - d)) + d_tube - a0;
            }else{
                y[i][0]= 0;
                y[i][JM-1] = d_tube;
            }
        }
        for (int j=1; j<JM-1;j++){ //create boundary nodes left/right
            x[0][j] = 0;
            y[0][j] = j*dy;
            x[IM-1][j] = l_tube;
            y[IM-1][j] = j*dy;
        }

        int stop = 0; int i = 0;
        while (stop == 0){      //identify idx_start (index at which neck starts to form)
            if (y[i][0] != y[0][0]) {
                idx_start = i-1;
                stop = 1;
            }
            i++;
        }

        double dy_temp;
        for (int i=1; i < IM-1 ; i++ ){ //initialize inner nodes
            dy_temp = (y[i][JM-1] - y[i][0] ) / (JM-1);
            for (int j=1; j<JM-1; j++){
                x[i][j] = x[i][0];
                y[i][j] = j*dy_temp + y[i][0];
            }
        }
        Print2FileXY("init");
    }//end of create_cosine_tube

    void SorensonComputeF(int method){
        ///Thomas & Middlecoff method
        const double dy = abs(y[0][JM-1] - y[0][0])/(JM-1);
        p_s = dy*p_s_percent;

        double f1_temp, f2_temp;
        double xxi, yxi, xxixi, yxixi, yeta, xeta;
        for (int i=1; i<IM-1; i++){         ///Calculate f1f2 in lower boundary wall
            xxi     = (1.0/2.0)*( x[i+1][0] - x[i-1][0] );      // d(x)/d(ksi)
            yxi     = (1.0/2.0)*( y[i+1][0] - y[i-1][0] );      // d(y)/d(ksi)
            xxixi   = x[i+1][0] + x[i-1][0] - 2.0*x[i][0];
            yxixi   = y[i+1][0] + y[i-1][0] - 2.0*y[i][0];

            //asumptions etc
            yeta  =  sqrt( pow(p_s*xxi,2) / (pow(xxi,2) + pow(yxi,2)) );  //its +-. We choose + so the Jacobean will be > 0  // d(y)/d(eta);
            xeta = - yxi*yeta/xxi;
            if (method == 2){
                yetaeta[i] = sorenson_relax*2.0*(y[i][1] - y[i][0]) + (1-sorenson_relax)*yetaeta[i];
                xetaeta[i] = sorenson_relax*2.0*(x[i][1] - x[i][0]) + (1-sorenson_relax)*xetaeta[i];
            }else if (method == 1){
                yetaeta[i] = 0;
                xetaeta[i] = 0;
            }else{
                std::cerr << " ERRROR IN SorensonComputeF, method number " << method << " is not recognised" << std::endl;
            }
            A[i][0] = pow(xeta ,2) + pow(yeta,2);
            C[i][0] = pow(xxi ,2) + pow(yxi ,2);
            J[i][0] = xxi*yeta - xeta*yxi;

            f1_temp = (1.0/pow(J[i][0], 3))*( - yeta*(A[i][0]*xxixi + C[i][0]*xetaeta[i]) + xeta*(A[i][0]*yxixi + C[i][0]*yetaeta[i]));
            f2_temp = (1.0/pow(J[i][0], 3))*( + yxi *(A[i][0]*xxixi + C[i][0]*xetaeta[i]) - xxi *(A[i][0]*yxixi + C[i][0]*yetaeta[i]));
            f_residual = (abs(f1_temp - f1[i][0]) + abs(f2_temp - f2[i][0]));

            f1[i][0] = f1_temp;
            f2[i][0] = f2_temp;

        }

        for (int i=1; i<IM-1; i++){ ///Project f1f2 in inner nodes up to top boundary wall(its symmetric to lower)
            for (int j=1; j<JM-1;j++){
                f1[i][j] = f1[i][0]*( exp(-p_tuner_f1*j) + exp(-p_tuner_f1*(JM-(j+1))) ) ;
                f2[i][j] = f2[i][0]*( exp(-p_tuner_f2*j) - exp(-p_tuner_f2*(JM-(j+1))) ) ; //if u use f[i][JM-1] change - to +
            }
        }

    }//end of SoresonComputeF

    double SurfaceTrap(){ ///na grapsw edw allo tropo evresis emvadou..
        std::vector<double> x_upper(IM, 0);
        std::vector<double> y_upper(IM, 0);
        for (int i=0; i<IM; i++){
            x_upper[i]= x[i][JM-1];
            y_upper[i]= y[i][JM-1] - d_tube/2;
        }
        return (2*trap(x_upper,y_upper) );
    }//end of SurfaceTrap

    void SurfaceJacobean(){ ///SOS we must have computed residual..so inner nodes have updated jacobean..

        surface = 0;
        for (int i=1; i<IM-1; i++){     //for inner nodes
            for (int j=1;j<JM-1; j++){
                surface = surface + J[i][j];
                if (J[i][j] < 0 ){
                    std::cout << "ERROR in AreaJacobean: on node   " << i <<  ", " << j << "  to prosimo tis J einai lathos " << std::endl;
                    kavalaei = 1;
                }

            }
        }

        double J_boundary ;
        for (int i=1; i<IM-1;i++){      //lower&upper boundary
            J_boundary = 0;
            J_boundary =    (1.0/2.0)*( x[i+1][0] - x[i-1][0] )*( y[i][1] - y[i][0] )  //1st order accuracy on eta
                            -(1.0/2.0)*( x[i][1] - x[i][0] )*( y[i+1][0] - y[i-1][0]);
            surface = surface + J_boundary/2;
            J[i][0] = J_boundary;
            if (J_boundary*J[2][2] < 0 ){
                std::cout << "ERROR in AreaJacobean: on boundary  " << i <<",0  J to prosimo tis J einai lathos " << std::endl;
                kavalaei = 1;
            }

            J_boundary = 0;
            J_boundary =    (1.0/2.0)*( x[i+1][JM-1] - x[i-1][JM-1] )*( y[i][JM-1] - y[i][JM-2] )  //1st order accuracy on eta
                            -(1.0/2.0)*( x[i][JM-1] - x[i][JM-2] )*( y[i+1][JM-1] - y[i-1][JM-1]);
            surface = surface + J_boundary/2;
            J[i][JM-1] = J_boundary;
            if (J_boundary*J[2][2] < 0 ){
                std::cout << "ERROR in AreaJacobean: on boundary  " << i <<",JM-1   to prosimo tis J einai lathos " << std::endl;
                kavalaei = 1;
            }
        }

        for (int j=1; j<JM-1;j++){      //left & right boundary
            J_boundary = 0;
            J_boundary =    (1.0/2.0)*( x[1][j] - x[0][j] )*( y[0][j+1] - y[0][j-1] )  //1st order accuracy on eta
                            -(1.0/2.0)*( x[0][j+1] - x[0][j-1] )*( y[1][j] - y[0][j]);
            surface = surface + J_boundary/2;
            J[0][j] = J_boundary;
            if (J_boundary*J[2][2] < 0 ){
                std::cout << "ERROR in AreaJacobean: on boundary 0, " << j <<" to prosimo tis J einai lathos " << std::endl;
                kavalaei = 1;
            }

            J_boundary = 0;
            J_boundary =    (1.0/2.0)*( x[IM-1][j] - x[IM-2][j] )*( y[IM-1][j+1] - y[IM-1][j-1] )  //1st order accuracy on eta
                            -(1.0/2.0)*( x[IM-1][j+1] - x[IM-1][j-1] )*( y[IM-1][j] - y[IM-2][j]);
            surface = surface + J_boundary/2;
            J[IM-1][j] = J_boundary;
            if (J_boundary*J[2][2] < 0 ){
                std::cout << "ERROR in AreaJacobean: on boundary IM-1, " << j <<" to prosimo tis J einai lathos " << std::endl;
                kavalaei = 1;
            }
        }

        ///corners
        J[0][0]         =    (1.0/1.0)*( x[1][0] - x[0][0] )*( y[0][1] - y[0][0] )  //1st order accuracy on eta
                            -(1.0/1.0)*( x[0][1] - x[0][0] )*( y[1][0] - y[0][0]);

        J[IM-1][0]      =    (1.0/1.0)*( x[IM-1][0] - x[IM-2][0] )*( y[IM-1][1] - y[IM-1][0] )  //1st order accuracy on eta
                            -(1.0/1.0)*( x[IM-1][1] - x[IM-1][0] )*( y[IM-1][0] - y[IM-2][0]);

        J[IM-1][JM-1]   =    (1.0/1.0)*( x[IM-1][JM-1] - x[IM-2][JM-1] )*( y[IM-1][JM-1] - y[IM-1][JM-2] )  //1st order accuracy on eta
                            -(1.0/1.0)*( x[IM-1][JM-1] - x[IM-1][JM-2] )*( y[IM-1][JM-1] - y[IM-2][JM-1]);

        J[0][JM-1]      =    (1.0/1.0)*( x[1][JM-1] - x[0][JM-1] )*( y[0][JM-1] - y[0][JM-2] )  //1st order accuracy on eta
                            -(1.0/1.0)*( x[0][JM-1] - x[0][JM-2] )*( y[1][JM-1] - y[0][JM-1 ]);
        if ( (J[0][0] < 0) || (J[IM-1][0] < 0) || (J[0][JM-1] < 0) || (J[IM-1][JM-1] < 0) ) {
                std::cout << "ERROR in AreaJacobean: a corner Boundary has wrong prosimo J " << std::endl;
                kavalaei = 1;
            }

        surface = surface + ( J[0][0] + J[IM-1][0] + J[IM-1][JM-1] + J[0][JM-1] )/4.0;
        Print2FileJ("final");
    }//end of SurfaceJacobean

    void UpdateMetrics(const int i, const int j){
        /// Calculate new metrics for specific node
        A[i][j] =    (1.0/4.0)*pow( (x[i][j+1] - x[i][j-1]) ,2)                     //g22
                    +(1.0/4.0)*pow( (y[i][j+1] - y[i][j-1]) ,2);
        B[i][j] =    (1.0/4.0)*(x[i+1][j] - x[i-1][j])*(x[i][j+1] - x[i][j-1])      //g12
                    +(1.0/4.0)*(y[i+1][j] - y[i-1][j])*(y[i][j+1] - y[i][j-1]);
        C[i][j] =    (1.0/4.0)*pow( (x[i+1][j] - x[i-1][j]) ,2)
                    +(1.0/4.0)*pow( (y[i+1][j] - y[i-1][j]) ,2);
        J[i][j] =    (1.0/4.0)*( x[i+1][j] - x[i-1][j] )*( y[i][j+1] - y[i][j-1] )  //g11
                    -(1.0/4.0)*( x[i][j+1] - x[i][j-1] )*( y[i+1][j] - y[i-1][j] );

    }//end of UpdateMetrics

    void UpdateResidual(double print_res){
        double res_x,res_x_cell;
        double res_y,res_y_cell;
        double xxi, xeta, xxixi, xetaeta, xxieta;
        double yxi, yeta, yxixi, yetaeta, yxieta;

        res_x = 0;
        res_y = 0;
        for (int j=1;j<JM-1;j++){
            for (int i=1;i<IM-1;i++){
                UpdateMetrics(i,j); //update metrics on specific node before creating the diagonal sys
                xxi     = (1.0/2.0)*( x[i+1][j] - x[i-1][j] );   // d(x)/d(ksi)
                xeta    = (1.0/2.0)*( x[i][j+1] - x[i][j-1] );
                xxixi   = x[i+1][j] + x[i-1][j] - 2.0*x[i][j];
                xetaeta = x[i][j+1] + x[i][j-1] - 2.0*x[i][j];
                xxieta  = (1.0)/(4.0)*( x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1] );
                res_x_cell = abs(
                                     A[i][j]*xxixi - 2.0*B[i][j]*xxieta + C[i][j]*xetaeta
                                    +pow(J[i][j], 2)*( f1[i][j]*xxi + f2[i][j]*xeta )
                                );
                yxi     = (1.0/2.0)*( y[i+1][j] - y[i-1][j] );   // d(y)/d(ksi)
                yeta    = (1.0/2.0)*( y[i][j+1] - y[i][j-1] );
                yxixi   = y[i+1][j] + y[i-1][j] - 2.0*y[i][j];
                yetaeta = y[i][j+1] + y[i][j-1] - 2.0*y[i][j];
                yxieta  = (1.0)/(4.0)*( y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1] );
                res_y_cell = abs(
                                     A[i][j]*yxixi - 2.0*B[i][j]*yxieta + C[i][j]*yetaeta
                                    + pow(J[i][j], 2)*( f1[i][j]*yxi + f2[i][j]*yeta )
                                );

                res_x = res_x + res_x_cell;
                res_y = res_y + res_y_cell;
                if (print_res == 1){
                    std::cout << "res [" << i << " ," << j << " ]"<< res_x_cell + res_y_cell << std::endl;
                }
            }
        }
        if (print_res == 1){
            std::cout << "res_x " << res_x << std::endl;
            std::cout << "res_y " << res_y << std::endl;
        }
        residual = res_x+res_y;
    }//end of UpdateResidual

    void HandleEndOfIteration(int iteration){
        UpdateResidual(0);
        //cout << sorenson_relax<< endl;
        for (int j = 0; j<JM; j++){      //update left/right boundaries to match with xi lines
                y[0][j] = y[idx_start/2][j];
                y[IM-1][j] = y[(IM-1) - idx_start/2][j];
            }
        if (p_tuner == TunerSorenson2){
            SorensonComputeF(TunerSorenson2);
            //update relaxation
            if (sorenson_relax < sorenson_relax_final){
                sorenson_relax = sorenson_relax + sorenson_multiplier*iteration*sorenson_relax;
                sorenson_counter++;
            }
            if (sorenson_relax > sorenson_relax_final){
                sorenson_relax = sorenson_relax_final;
            }
        }
    }

    void PoisonGridPointJacobi(const double res_tol, const double max_iterations,
                            const double relaxation, const int info_counter){
        int iteration =0; int counter = 0; //to not print every iteration
        double x_new, y_new;

        while ( ( (iteration < max_iterations) &&  (residual > res_tol ) ) || (sorenson_relax < sorenson_relax_final) )
        {
            counter++;      iteration++;
            if (counter == info_counter){
                counter = 0;
                std::cout <<"iteration: " << iteration << "   Residual: " << residual << std::endl;
                if (p_tuner == TunerSorenson2){
                                std::cout << "f_residual " << f_residual << std::endl;
                }
            }


            for (int i=1; i<IM-1; i++){
                for (int j=1; j<JM-1; j++){
                    UpdateMetrics(i,j);
                    x_new =(
                                -A[i][j]*(x[i+1][j] + x[i-1][j])
                                -C[i][j]*(x[i][j+1] + x[i][j-1])
                                +0.5*B[i][j]*( x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1] )
                                -0.5*pow(J[i][j], 2)*f1[i][j]*( x[i+1][j] - x[i-1][j])
                                -0.5*pow(J[i][j], 2)*f2[i][j]*( x[i][j+1] - x[i][j-1])
                            )/(-2.0*(A[i][j] + C[i][j]));

                    y_new =(
                                -A[i][j]*(y[i+1][j] + y[i-1][j])
                                -C[i][j]*(y[i][j+1] + y[i][j-1])
                                +0.5*B[i][j]*( y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1] )
                                -0.5*pow(J[i][j], 2)*f1[i][j]*( y[i+1][j] - y[i-1][j])
                                -0.5*pow(J[i][j], 2)*f2[i][j]*( y[i][j+1] - y[i][j-1])
                            )/(-2.0*(A[i][j] + C[i][j]));

                    x[i][j] = relaxation*x_new + (1-relaxation)*x[i][j];
                    y[i][j] = relaxation*y_new + (1-relaxation)*y[i][j];

                }
            }
            HandleEndOfIteration(iteration);
        }

        std::cout << std::endl;
        if (residual < res_tol){
            std::cout << "POISSON POINT JACOBI SOLVER: Computation of Poisson grid finished successfully. Residual = " << std::endl;
            std::cout << "->Residual = " << residual <<std:: endl;
            if (p_tuner == TunerSorenson2){
             //   std::cout << "->F_residual = " << residual << " ...f_relaxation after " << f_relax_it << " th iteration: " << f_relax <<std:: endl;
            }
            std::cout << "->Iteration = " << iteration << std::endl;
        }else{
            std::cout << "POISSON POINT JACOBI SOLVER: Computation of Poisson grid was NOT successful. Max iteration limit reached ( " << max_iterations << " )" << std::endl;
            std::cout << "->Residual = " << residual << std::endl;
        }
        Print2FileXY("final");
    }// end of Point Jacobi

    void PoisonGridImplicit(const double res_tol, const double max_iterations,
                            const double relaxation, const int info_counter){

        int counter = 0; //to not print every iteration
        int iteration = 0;
        while ( ( (iteration < max_iterations) &&  (residual > res_tol ) )   || (sorenson_relax < sorenson_relax_final) )
        {
            counter++;      iteration++;
            if (counter == info_counter){
                counter = 0;
                std::cout <<"iteration: " << iteration << "   Residual: " << residual << std::endl;
                if (p_tuner == TunerSorenson2){
                                std::cout << "f_relax " << sorenson_relax << std::endl;
                }
            }

            /// NORTH - SOUTH sweep
            std::vector<double> b(IM,0), RHSx(IM,0), RHSy(IM,0), x_new(IM,0), y_new(IM,0); //because we will do only x sweeps
            std::vector<double> a(IM-1,0),c(IM-1,0);
            for (int j=1;j<JM-1;j++){
                for (int i=1;i<IM-1;i++){
                    //Create diagonal sys
                    UpdateMetrics(i,j); //update metrics on specific node before creating the diagonal sys
                    b[i]    =   -2.0*(A[i][j] + C[i][j]);
                    a[i-1]  =   A[i][j] - 0.5*pow(J[i][j], 2)*f1[i][j];   //this changes from laplace to poisson
                    c[i]    =   A[i][j] + 0.5*pow(J[i][j], 2)*f1[i][j];   //this changed from laplace to poisson
                    RHSx[i] =    (0.5)*B[i][j]
                                *( x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1] )
                                - C[i][j] *(x[i][j+1] + x[i][j-1])
                                - (pow(J[i][j], 2) * f2[i][j]* 0.5*(x[i][j+1] - x[i][j-1]) );
                    RHSy[i] =   +(0.5)*B[i][j]
                                *( y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1] )
                                - C[i][j] *(y[i][j+1] + y[i][j-1])
                                - (pow(J[i][j], 2) * f2[i][j]* 0.5*(y[i][j+1] - y[i][j-1]) );

                }
                //enforce boundaries(i=0, i=IM-1
                b[0]=1; b[IM-1]=1; c[0]=0; a[IM-2]=0;
                RHSx[0]=x[0][j];   RHSx[IM-1]=x[IM-1][j];
                RHSy[0]=y[0][j];   RHSy[IM-1]=y[IM-1][j];
                //solve row
                tridiagonal_solver(a,b,c,RHSx,x_new);   //Solve for x
                tridiagonal_solver(a,b,c,RHSy,y_new);   //Solve for x
                //Update row
                for (int i=1; i<IM-1; i++){
                    x[i][j] = x_new[i]*relaxation + (1-relaxation)*x[i][j];
                    y[i][j] = y_new[i]*relaxation + (1-relaxation)*y[i][j];
                }
            }
            ///WEST - EAST sweep
            b.resize(JM), RHSx.resize(JM), RHSy.resize(JM), x_new.resize(JM), y_new.resize(JM);
            a.resize(JM-1), c.resize(JM-1);

            for (int i=1;i<IM-1;i++){
                for (int j=1;j<JM-1;j++){
                    //Create diagonal sys
                    UpdateMetrics(i,j); //update metrics on specific node before creating the diagonal sys
                    b[j]    =   -2.0*(A[i][j] + C[i][j]);
                    a[j-1]  =   C[i][j] - 0.5*pow(J[i][j], 2)*f2[i][j];   //this changes from laplace to poisson
                    c[j]    =   C[i][j] + 0.5*pow(J[i][j], 2)*f2[i][j];   //this changed from laplace to poisson
                    RHSx[j] =    (0.5)*B[i][j]
                                *( x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1] )
                                - A[i][j] *(x[i+1][j] + x[i-1][j])
                                - (pow(J[i][j], 2) * f1[i][j]* 0.5*(x[i+1][j] - x[i-1][j]) );
                    RHSy[j] =    (0.5)*B[i][j]
                                *( y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1] )
                                - A[i][j] *(y[i+1][j] + y[i-1][j])
                                - (pow(J[i][j], 2) * f1[i][j]* 0.5*(y[i+1][j] - y[i-1][j]) );
                }
                //enforce boundaries(j=0, j=JM-1)
                b[0]=1; b[JM-1]=1; c[0]=0; a[JM-2]=0;
                RHSx[0]=x[i][0];   RHSx[JM-1]=x[i][JM-1];
                RHSy[0]=y[i][0];   RHSy[JM-1]=y[i][JM-1];
                //solve row
                tridiagonal_solver(a,b,c,RHSx,x_new);   //Solve for x
                tridiagonal_solver(a,b,c,RHSy,y_new);   //Solve for x
                //Update row
                for (int j=1; j<JM-1; j++){
                    x[i][j] = x_new[j]*relaxation + (1-relaxation)*x[i][j];
                    y[i][j] = y_new[j]*relaxation + (1-relaxation)*y[i][j];
                }
            }
            HandleEndOfIteration(iteration);
        }
        std::cout << std::endl;
        if (residual < res_tol){
            std::cout << "POISSON IMPLICIT SOLVER: Computation of Poisson grid finished successfully." << std::endl;
            std::cout << "->Residual = " << residual <<std:: endl;
            if (p_tuner == TunerSorenson2){
               // std::cout << "->F_residual = " << residual << " ...f_relaxation after " << f_relax_it << " th iteration: " << f_relax <<std:: endl;
            }
            std::cout << "->Iteration = " << iteration << std::endl;
        }else{
            std::cout << "POISSON IMPLICIT SOLVER: Computation of Poisson grid was NOT successful. Max iteration limit reached ( " << max_iterations << " )" << std::endl;
            std::cout << "->Residual = " << residual << std::endl;
        }
        Print2FileXY("final");
    }//end of PoissonGridImplicit
};

///
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///

int main()
{
    /// USER INPUTS
    //Cosine tube characteristics
    const double L      =  5;      const double D      = 4;
    const double l_neck =  4;      const double d_neck = 1;
    const double x_nodes_num =  101;  const double y_nodes_num = 121;
    GRid2DCosineTube grid(x_nodes_num,y_nodes_num, L, D, l_neck, d_neck);
    //Characteristics of tuning
    const double tuner = grid.TunerSorenson1; //1 for thomas & middlecoff, 2 for sorenson
    const double ds_percentage=0.6; // percentage of dy
    const double tuner_f1 = 0.1; const double tuner_f2=0.1;
    const double frelax_s = 0.01; const double frelax_f = 0.3;
    const double sorenson_multiplier = 0.001; //slows the sorenson relaxation update

    //Solver Characteristics
    const double residual_tol = pow(10, -4); const double relax = 1;
    const int max_it = 1000; //max iterations
    const int info_counter=100; // how many iterations do u want to pass when u get information about the residual
    const int solver = grid.SolverImplicit;
    grid.ComputeF(tuner, ds_percentage,tuner_f1, tuner_f2, frelax_s, frelax_f, sorenson_multiplier); //update f1 f2


    /// Solve
    grid.solve(solver, residual_tol, max_it, relax, info_counter);

    ///Post Process
    grid.PostProcess();


    std::cout << "END OF SCRIPT"  << std::endl;
    return 0;
}
