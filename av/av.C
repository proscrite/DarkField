/**
 * avalanche.cc
 *
*/
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <TCanvas.h>
#include <TApplication.h>
#include <TFile.h>

#include "MediumMagboltz.hh"
#include "ComponentElmer.hh"
#include "Sensor.hh"
#include "ViewField.hh"
#include "Plotting.hh"
#include "ViewFEMesh.hh"
#include "ViewSignal.hh"
#include "GarfieldConstants.hh"
#include "Random.hh"
#include "AvalancheMicroscopic.hh"

using namespace Garfield;

int ie_count = 0;
void hd_inelastic(double x, double y, double z, double t, int type, int level, Medium* m) {

    // Output the collision information: note type 4 is "excitation"
    std::cout << x << " " << y << " " << z << " " << t << " " << type << " " << level << std::endl;
    ie_count += 1;
}

int st_count = 0;
void hd_step(double x, double y, double z, double t, double e, double dx, double dy, double dz, bool hole) {

    // Output the step information
    if(st_count % 1000 == 0) {
      std::cout << x << " " << y << " " << z << " " << t << " " << -1 << " " << -1 << std::endl;
    }
    st_count += 1;
}

int main(int argc, char * argv[]) {

    TApplication app("app", &argc, argv);

    const int nelectrons = 100.;

    // Set relevant parameters.
    //const double d_anode1 = 5;
    //const double d_anode2 = 5;
    const double d_wall = 1.27;
    //const double d_hlwire = 1;
    const double d_rwire = 0.0125;

    const double axis_x = 0.5;
    const double axis_y = d_wall;
    const double axis_z = 0.5; //fmax(d_anode1,d_anode2);

    // ---------------------------------------------------------------------------------------------------------------
    // Create several canvases for the plots.
    TCanvas * cGeom = new TCanvas("geom","Geometry/Avalanche/Fields");
    //TCanvas * cSignal = new TCanvas("signal","Signal");

    // Define the medium.
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetTemperature(293.15);                  // Set the temperature (K)
    gas->SetPressure(7500.);                       // Set the pressure (Torr)
    gas->EnableDrift();                           // Allow for drifting in this medium
    gas->SetComposition("xe", 100.);   // Specify the gas mixture

    // Import an Elmer-created LEM and the weighting field for the readout electrode.
    ComponentElmer * elm = new ComponentElmer("wire/mesh.header","wire/mesh.elements","wire/mesh.nodes",
    "wire/dielectrics.dat","wire/wire.result","cm");
    elm->EnablePeriodicityX();
    elm->EnablePeriodicityY();
    elm->SetMedium(0,gas);

    // Set up a sensor object.
    Sensor* sensor = new Sensor();
    sensor->AddComponent(elm);
    sensor->SetArea(-1*axis_x,-1*axis_y,-1*axis_z,axis_x,axis_y,axis_z);

    // Create an avalanche object
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    aval->SetCollisionSteps(100);
    aval->EnableSignalCalculation();
    aval->SetUserHandleInelastic(hd_inelastic);
    aval->SetUserHandleStep(hd_step);

    // Set up the object for field visualization.
    ViewField * vf = new ViewField();
    vf->SetSensor(sensor);
    vf->SetCanvas(cGeom);
    vf->SetArea(-1*axis_x,-1*axis_y,axis_x,axis_y);
    vf->SetNumberOfContours(80);
    vf->SetNumberOfSamples2d(50,50);
    vf->SetPlane(0,-1,0,0,0,0);

    // Set up the object for drift line visualization.
    ViewDrift* viewDrift = new ViewDrift();
    viewDrift->SetArea(-1*axis_x,-1*axis_y,-1*axis_z,axis_x,axis_y,axis_z);
    aval->EnablePlotting(viewDrift);

    // Set up the object for FE mesh visualization.
    ViewFEMesh * vFE = new ViewFEMesh();
    vFE->SetCanvas(cGeom);
    vFE->SetComponent(elm);
    vFE->SetPlane(-1,0,0,0,0,0);
    vFE->SetFillMesh(true);
    vFE->SetColor(1,kGray);
    //vFE->SetColor(2,kYellow+3);
    //vFE->SetColor(3,kYellow+3);

    // Calculate the avalanche.
    for (int i = 0; i < nelectrons; i++) {

        // Set the electron start parameters.
        double ri = d_rwire+1e-6;
        double thetai = -RndmUniform()*TwoPi/2;
        double xi = 0;
        double yi = ri*cos(thetai);
        double zi = ri*sin(thetai);

        // For obtaining the endpoint:
        double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1;
        int status;

        ie_count = 0;
        st_count = 0;
        std::cout << "Beginning avalanche " << i << " of " << nelectrons << " from ("
                  << xi << ", " << yi << ", " << zi << ")" << std::endl;
        std::cout << "# x y z t type level" << std::endl;
        aval->AvalancheElectron(xi, yi, zi, 0., 0., 0., 0., 0.);
        std::cout << "... avalanche complete with " << ie_count
                  << " inelastics and " << aval->GetNumberOfElectronEndpoints()
                  << " electron tracks." << std::endl;
        aval->GetElectronEndpoint(0,x0,y0,z0,t0,e0,x1,y1,z1,t1,e1,status);
        std::cout << " -- Electron start (x,y,z,t,e): ("
                  << x0 << ", " << y0 << ", " << z0 << ", " << t0 << ", " << e0
                  << ") and end: ("
                  << x1 << ", " << y1 << ", " << z1 << ", " << t1 << ", " << e1
                  << ")" << std::endl;

    }

    // ---------------------------------------------------------------------------------------------------------------

    /*// Create ROOT histograms of the signal and a file in which to store them.
    TFile * f = new TFile("avalanche_signals.root","RECREATE");
    TH1F * hS = new TH1F("hh","hh",nsBins,0,tEnd);               // total signal
    TH1F * hInt = new TH1F("hInt","hInt",nsBins,0,tEnd);         // integrated signal

    // Fill the histograms with the signals.
    //  Note that the signals will be in C/(ns*binWidth), and we will divide by e to give a signal in e/(ns*binWidth).
    //  The total signal is then the integral over all bins multiplied by the bin width in ns.
    for(int i = 0; i < nsBins; i++) {
      double wt = sensor->GetSignal("wtlel",i)/ElementaryCharge;
      sum += wt;
      hInt->Fill(i*bscale,sum);
      hS->Fill(i*bscale,wt);
    }

    // Write the histograms to the TFile.
    hS->Write();
    hInt->Write();
    f->Close();*/

    // Fields
    double z0 = 2.5;
    double ex, ey, ez, vv;
    Medium * tmed;
    int tstatus = 0;
//    std::cout << "# Electric field x  y  z  Ex  Ey  Ez  V status" << std::endl;
//    for(double zz = -1*z0*2; zz < z0*2; zz += z0/1000) {
//        double xy = 0.0;
//        sensor->ElectricField(xy,xy,zz,ex,ey,ez,vv,tmed,tstatus);
        //std::cout << xy << " " << xy << " " << zz << " " << ex << " " << ey
        //          << " " << ez << " " << vv << " " << tstatus << std::endl;
//    }

    // ---------------------------------------------------------------------------------------------------------------
    // Create plots.
    vFE->SetArea(-1*axis_y,-1*axis_z,0,axis_y,axis_z,0);  // note: here the x-y axes correspond to projection chosen
                                                            //       z-axis is irrelevant for 2D projections
    //vf->PlotContour("v"); // uncomment this to plot the contours of the potential
    //vSignal->PlotSignal("wtlel");

    vFE->EnableAxes();             // comment this to disable creation of independent axes when contours are plotted
    vFE->SetViewDrift(viewDrift);  // comment this to remove the avalanche drift line from the plot when contours are plotted
    vFE->SetXaxisTitle("z (cm)");
    vFE->SetYaxisTitle("y (cm)");
    vFE->Plot();

    app.Run(kTRUE);

    return 0;
}
