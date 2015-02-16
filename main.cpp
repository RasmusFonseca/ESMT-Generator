#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>

using namespace std;

int printUsage(char* progName)
{
    cout<<"Usage: "<<progName<<" <pointset type> [options]"<<endl<<endl;

    cout<<"Pointset type must be one of:"<<endl;
    cout<<"\tcube : points distributed randomly in a unit cube"<<endl;
    cout<<"\tsphere : points distributed randomly in a unit sphere"<<endl;
    cout<<"\tedge : points distributed randomly on diagonal of a hypercube"<<endl;
    cout<<"\tspheresurface : points distributed evenly on a unit sphere surface"<<endl;
    cout<<"\tsausage : the corners of a sequence of regular simplices that are face-to-face"<<endl;
    cout<<"\tgrid : the vertices of a cubic grid of width floor(n^(1/d))"<<endl;
    cout<<endl;
    cout<<"Options can be any list of:"<<endl;
    cout<<"\t-n <int> : number of points (standard=10)"<<endl;
    cout<<"\t-d <int> : dimension of points (standard=2)"<<endl;
    cout<<"\t-s <int> : seed for generating random points (standard=time())"<<endl;
    cout<<"\t-name <string> : name of point set (standard=\"\")"<<endl;

    return -1;
}

void printIntro(int n, int d, string &name);

void printCube(int n, int d, string &name);
void printEdge(int n, int d, string &name);
void printSphere(int n, int d, string &name);
void printSphereSurface(int n, int d, string &name);
void printSausage(int n, int d, string &name);
void printGrid(int n, int d, string &name);

void printOutro();

int main(int argc, char** argv)
{
    if(argc<=1) return printUsage(argv[0]);

    int n = 1000;
    int d = 3;
    unsigned int seed = (unsigned int)time(NULL);
    string name = "";
    for(int i=2;i<argc;i++){
        if(string(argv[i])=="-n") { n = atoi(argv[++i]); continue; }
        if(string(argv[i])=="-d") { d = atoi(argv[++i]); continue; }
        if(string(argv[i])=="-s") { seed = atoi(argv[++i]); continue; }
        if(string(argv[i])=="-name") { name = argv[++i]; continue; }
    }
    srand(seed);

    if(string(argv[1])=="cube") {           printCube(n, d, name);      		return 0; }
    if(string(argv[1])=="sphere") {         printSphere(n, d, name);    		return 0; }
    if(string(argv[1])=="spheresurface") {  printSphereSurface(n, d, name); 	return 0; }
    if(string(argv[1])=="sausage") {        printSausage(n, d, name);   		return 0; }
    if(string(argv[1])=="edge") {        	printEdge(n, d, name);   			return 0; }
    if(string(argv[1])=="grid") {           printGrid(n, d, name);   			return 0; }

    cerr<<"Unknown point set type: "<<argv[1]<<endl;
    return -1;
}

void printIntro(int n, int d, string &name)
{
    cout<<"33D32945 STP File, STP Format Version 1.0"<<endl;
    cout<<""<<endl;
    cout<<"SECTION Comments"<<endl;
    cout<<"Name \""<<name<<"\""<<endl;
    cout<<"Problem \"Euclidean Steiner Tree Problem\""<<endl;
    cout<<"Remark \"Generated with STPGenerator. "<<d<<" dimensions.\""<<endl;
    cout<<"END"<<endl;
    cout<<""<<endl;
    cout<<"SECTION Graph"<<endl;
    cout<<"Nodes "<<n<<endl;
    cout<<"END"<<endl;
    cout<<""<<endl;
    cout<<"SECTION Coordinates"<<endl;
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double fNorm(double mu, double sigma)
{
    //See http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    double rand1 = rand() / ((double) RAND_MAX);
    if(rand1 < 1e-100) rand1 = 1e-100;
    rand1 = -2 * log(rand1);
    double rand2 = (rand() / ((double) RAND_MAX)) * 2*M_PI;

    return (sigma * sqrt(rand1) * cos(rand2)) + mu;

}

void printCube(int n, int d, string &name)
{
    printIntro(n, d, name);
    for(int i=0;i<n;i++){
        for(int j=0;j<d;j++){
            cout<<"D";
        }
        for(int j=0;j<d;j++){
            printf(" %.10f",fRand(0.0,1.0));
        }
        cout<<endl;
    }
    printOutro();
}

double lenSq(vector<double> &vec){
    double ret = 0.0;
    for(size_t i=0;i<vec.size();i++){
        ret+=vec[i]*vec[i];
    }
    return ret;
}

void printSphere(int n, int d, string &name)
{
    printIntro(n, d, name);
    for(int i=0;i<n;i++){
        for(int j=0;j<d;j++){
            cout<<"D";
        }
        vector<double> coords(d,1.0);
        while(lenSq(coords)>1.0){
            for(int j=0;j<d;j++){
                coords[j] = fRand(-1.0,1.0);
            }
        }
        for(int j=0;j<d;j++){
            printf(" %.10f",coords[j]);
        }
        cout<<endl;
    }
    printOutro();

}

/** Randomly distributed points on sphere surface */
void printSphereSurface(int n, int d, string &name)
{

    printIntro(n, d, name);
    for(int i=0;i<n;i++){
        for(int j=0;j<d;j++){
            cout<<"D";
        }
        vector<double> coords(d,1.0);
        double len;
        while( (len=lenSq(coords))>1.0){
            for(int j=0;j<d;j++){
                coords[j] = fRand(-1.0,1.0);
            }
        }
        for(int j=0;j<d;j++){
            coords[j] /= sqrt(len);
        }

        for(int j=0;j<d;j++){
            printf(" %.10f",coords[j]);
        }
        cout<<endl;
    }
    printOutro();
}

void printEdge(int n, int d, string &name)
{
    printIntro(n, d, name);
    for(int i=0;i<n;i++){
        for(int j=0;j<d;j++){
            cout<<"D";
        }
        vector<double> coords(d,0.0);
        double coord = fRand(0.0,1.0);
        for(int j=0;j<d;j++){
            coords[j] = coord;
        }

        for(int j=0;j<d;j++){
            printf(" %.10f",coords[j]);
        }
        cout<<endl;
    }
    printOutro();
}

void printClusters(int n, int d, string &name)
{
    printIntro(n, d, name);

    vector<double> clusterCenter(d,0.0);

    for(int i=0;i<n;i++){
        if(i%(n/4)==0)
            for(int j=0;j<d;j++) clusterCenter[j] = fRand(0.0, 1.0);


        for(int j=0;j<d;j++){
            cout<<"D";
        }

        vector<double> coords(d,0.0);
        for(int j=0;j<d;j++){
            coords[j] = fNorm(clusterCenter[j], 0.05);
        }

        for(int j=0;j<d;j++){
            printf(" %.10f",coords[j]);
        }
        cout<<endl;
    }
    printOutro();
}



void printSausage(int n, int d, string &name)
{
    printIntro(n, d, name);

    int corners = d+1;
    double angle = -1.0/d;
    vector< vector<double> > prevPoints(d+1, vector<double>(d, 0.0));

    //Initialize prevPoints with a regular d-dimensional simplex
    double tmp = 0.0;
    for(int i=0;i<d;i++){
        double nCoord = sqrt(1-tmp);
        prevPoints[i][i] = nCoord;
        double rCoord = (angle - tmp)/nCoord;

        for(int s=i+1;s<corners;s++)
            prevPoints[s][i] = rCoord;

        tmp += rCoord*rCoord;
    }

    //Extend the simplex by consecutively pushing apex point through base
    for(int i=d+1;i<n;i++){
        vector<double> newPoint(d, 0.0);
        for(int j=0;j<d;j++){
            double center_j = (prevPoints[i-1][j] + prevPoints[i-2][j] + prevPoints[i-3][j])/3.0;
            newPoint[j] = center_j + (center_j-prevPoints[i-4][j]);
        }
        prevPoints.push_back(newPoint);

    }

    //Print n points of sausage
    for(int i=0;i<n;i++){

        for(int j=0;j<d;j++){
            cout<<"D";
        }

        for(int j=0;j<d;j++){
            printf(" %.10f",prevPoints[i][j]);
        }
        cout<<endl;
    }

    printOutro();
}


void printGrid(int n, int d, string &name)
{
    printIntro(n, d, name);
    int width = (int)(pow(n,1.0/d)+0.00001);
    vector<int> gridPos(d,0);
    for(int i=0;i<pow(width,d);i++){
        for(int j=0;j<d;j++){
            cout<<"D";
        }
        for(int j=0;j<d;j++){
            printf(" %.1f",gridPos[j]*1.0);
        }
        gridPos[0]++;
        for(int j=0;j<d-1;j++)
            if(gridPos[j]==width){
                gridPos[j+1]++;
                gridPos[j]=0;
            }
        cout<<endl;
    }
    printOutro();
}

void printOutro()
{
    cout<<"END"<<endl;
}
