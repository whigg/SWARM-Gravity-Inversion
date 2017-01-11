#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


#include "coord.h"
#include "numrs.h"


using namespace std;


//#define MAXLINE 300

//short int ACCURACY = 1;

//int LAG_ORD = 8;

//int NEOP;
//double *EOPMT;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int main (int argc, char *argv[])
{
    int year, month, day, hour, min, mjds, mjde,  gps_i, i, n, GPS_S;
    double JD0, jdi, sec, xe[3], ve[3],  xi[3], vi[3] = {0},  tt, 
        tjd[2], viv[3], vix[3];
        
//    char line[MAXLINE];  
    string stdname, card, fsp3, fgpsi, fgpse, feop;  
    istringstream iline;

    InfStruct info;

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input file name>" << endl;
        return EXIT_FAILURE;
    }
    ifstream if_std, if_sp3;
    ofstream of_pvi, of_pve; 
    if_std.open(argv[1]); 
//    ostringstream line;
//    stringstream sl   ine; 
    string line;
    while(!if_std.eof()) {
        getline(if_std, line);
//        cout << line << endl; 
        iline.clear(); iline.str(line);
        iline >> card;
        if (card == "INPUT") {iline >> fsp3; if_sp3.open(fsp3.c_str());}
        if (card == "OUTPUT") {iline >> fgpsi >> fgpse; of_pvi.open(fgpsi.c_str()); of_pve.open(fgpse.c_str()); }
        if (card == "EOP") {iline >> feop;}
        if (card == "EPOCH") {iline >> year >> month >> day; }
    }
    
    JD0 = julian_date ((short)year,(short)month,(short)day,0);
    GPS_S = (int)((JD0 - T0) * 86400 + 0.5);

    mjds = (int)(JD0 - 2400000.5);
    mjde = mjds + 1;
    eop_open (feop.c_str(), mjds - 1, mjde + 1);
    

    if_sp3.ignore(2);
    char flag;
    if_sp3.get(flag); // cout << flag;
    for (i = 0; i < 22; i++) {
        getline(if_sp3, line);
//        cout << line << endl;
    }

    of_pvi.precision(14);
    of_pve.precision(14);

    char label[4];
//    cout << flag << '\t' << JD0 << '\t'<< day << endl;
    while(!if_sp3.eof()) {
        if_sp3.get(label,4); 
        if (!strcmp (label, "EOF")) {
//            cout << label << endl;
            break;
        }
        getline(if_sp3, line);
//        cout << line << endl;
        iline.clear(); iline.str(line);
//        iline.ignore(3);
//        iline.get(label,4); 
//        cout << label << endl;
//        if (!strcmp (label, "EOF")) {cout << label << endl;break;}
//        if ( label == "EOF") {cout << label << endl;break;}
        iline >> year >> month >> day >> hour >> min >> sec;
        jdi = julian_date ((short)year,(short)month,(short)day,0);
        gps_i = (int)((jdi - T0) * 86400 + 0.5) + hour * 3600 + min * 60 + sec;

        getline(if_sp3, line);
//        cout << line;
        iline.clear(); iline.str(line);
        iline.ignore(4);
        iline >> xe[0] >> xe[1] >> xe[2];
        for (n = 0; n < 3; n ++) xe[n] = xe[n] * 1000;
    
        if (flag == 'V') {
            getline(if_sp3, line);
//            cout << line;
            iline.clear(); iline.str(line);
            iline.ignore(4);
            iline >> ve[0] >> ve[1] >> ve[2];
    
            for (n = 0; n < 3; n ++) ve[n] = ve[n] * 0.1;
            
        }


        tt = gps_i - GPS_S + 19 + 32.184;
        tjd[0] = JD0;    tjd[1] = tt / 86400.0;
        getinfo (tjd, 2, &info);

        brmul(info.c_ei, xe, 3,3,1, xi);

        if (flag == 'V') {
            brmul(info.c_ei, ve, 3,3,1, viv);
            brmul(info.c_eidot, xe, 3,3,1, vix);
            for (n = 0; n < 3; n ++) vi[n] = viv[n] + vix[n];
        }

        of_pve << gps_i << '\t' << xe[0] << '\t' << xe[1] << '\t' << xe[2] << '\t';
        of_pvi << gps_i << '\t' << xi[0] << '\t' << xi[1] << '\t' << xi[2] << '\t';
        if (flag == 'V') {
            of_pve << ve[0] << '\t' << ve[1] << '\t' << ve[2];
            of_pvi << vi[0] << '\t' << vi[1] << '\t' << vi[2];
        }
        of_pve << endl;
        of_pvi << endl;

    }
    
    if_std.close();
    if_sp3.close();
    of_pve.close();
    of_pvi.close();


    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 
