typedef struct worker_struct_in
{
    CMPLX *p_srcAs;     
    double *p_bessel_table; 
    int M, N, nCha, nPHS, r_range, loopInd;
    float delta_r, delta_rho;
    worker_struct_in(CMPLX *a=NULL, double *b=NULL, int c=0, int d=0, int e=0, int f=0, int g=0, float h=0, float i=0, int ind=0): p_srcAs(a), p_bessel_table(b), M(c), N(d), nCha(e), nPHS(f), r_range(g), delta_r(h), delta_rho(i), loopInd(ind){}
} worker_struct_in;


typedef struct worker_struct_out
{
    CMPLX *p_dstAs;
    worker_struct_out(CMPLX *a=NULL): p_dstAs(a){}
} worker_struct_out;


bool calc_bessel_table(std::vector<double> &bessel_table, int M, int N, int nPHS, float delta_r, float delta_rho, int r_range)
{
	//parameters for Hankel Transform
	//M: Number of samples on a spoke
	//N: total number of orders
	//rho: a 1D vector contains rho values (in source image)
	//r:   a 1D vector contains r values   (for final image)
	//r_range: length of vector r
	int	   z_ord = (int) N/2;
	float  pii  = 3.14159265f;	// f means float. to avoid error

    stringstream ss;
    ss << getenv("SimMeasData") << "/bessel_table/"<< N << "_" << M << "_" << r_range << ".bin";
    std::string filename = ss.str();
    ICE_OUT("bessel table filename: " << filename);
    std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
    if( f.good() && f.is_open())
    {
         f.read((char*)&bessel_table[0], sizeof(double)*N*M*r_range);
         f.close();
    }
    else
    {
        std::cout << "bessel table file not found, creating it..." << std::endl;
        std::ofstream fo(filename.c_str(), std::ios::out | std::ios::binary);
        if (fo.is_open() == false)
        {
            ICE_OUT("cannot open file " << filename.c_str() );
            return false;
        }

        for (int phase = 0; phase < nPHS; phase++)
        {	
            for(int ord = -(N/2); ord < (N/2); ord++)		//lines
            {	
                for(int rr = 0; rr < r_range; rr++)	//columns
                {
                    float r=(float)(rr+1)*delta_r;
                    for(int rhoo = 0; rhoo < M; rhoo++)
                    {
                        float rho = delta_rho/2 + (float)rhoo*delta_rho;                            
                        bessel_table[(ord+z_ord)*M*r_range  + rhoo*r_range + rr] = bessj(ord , (float)(pii*rho*r));
                    }
                }
            }
        }
        fo.write((char*)&bessel_table[0], sizeof(double)*N*M*r_range);
        fo.close(); 
    }
    return true;
}



bool Hankel_Transform_Parallel(worker_struct_in in , worker_struct_out out)
{
    CMPLX *p_dstAs = out.p_dstAs;
    CMPLX *p_srcAs = in.p_srcAs;
    double *p_bessel_table = in.p_bessel_table;
    int M = in.M;
    int N = in.N;
    int nCha = in.nCha;
    int nPHS = in.nPHS;
    int r_range = in.r_range;
    float delta_r = in.delta_r;
    float delta_rho = in.delta_rho;
    
    int loopInd = in.loopInd;
    
	//parameters for Hankel Transform
	//M: Number of samples on a spoke
	//N: total number of orders
	//rho: a 1D vector contains rho values (in source image)
	//r:   a 1D vector contains r values   (for final image)
	//r_range: length of vector r
	int	   z_ord = (int) N/2;
 

	// Hankel Formula //
	float my_j;
	float r;
	float rho;
	complex<float> temp (0.0 , 0.0);
	complex<float> fn_r (0.0 , 0.0);
	complex<float> immm(0, 1);
	int ipow ;

    p_srcAs += loopInd * M;
    p_dstAs += loopInd * r_range;
	for (int phase = 0; phase < nPHS; phase++)
	{	
        for(int ord = -(N/2); ord < (N/2); ord++)	//lines
		{
            for (int rr = 0; rr < r_range; rr++) // columns
            {
                r = (float)(rr + 1) * delta_r;
                fn_r = complex<float>(0.0, 0.0);

                for (int rhoo = 0; rhoo < M; rhoo++)
                {
                    rho = delta_rho / 2 + (float)rhoo * delta_rho;
                    my_j = p_bessel_table[(ord + z_ord) * M * r_range + rhoo * r_range + rr];

                    temp = complex<float>(rho * my_j * p_srcAs->re, rho * my_j * p_srcAs->im);
                    fn_r += temp;
                    p_srcAs++;
                }
                // Pi gain (to increase contrast of final image)
                fn_r *= 2 * PI;

                // instead of (fn_r *= std::pow(immm,ipow) or exp(ipow*log(immm))) just because the (pow() ) doesn't work here :(
                //--------
                ipow = (ord + z_ord) % 4;
                if (ipow == 1)
                {
                    fn_r *= immm;
                }
                else if (ipow == 2)
                {
                    fn_r *= -1.0f;
                }
                else if (ipow == 3)
                {
                    fn_r *= -immm;
                }
                //--------------

                p_dstAs->re = real(fn_r);
                p_dstAs->im = imag(fn_r);
                p_dstAs++;
                if (rr != (r_range - 1))
                    p_srcAs -= M;
            }
            p_srcAs += (nCha - 1) * M;
            p_dstAs += (nCha - 1) * r_range;
        }
	}
	return true;
}