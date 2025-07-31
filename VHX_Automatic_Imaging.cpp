#include <windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "include/VHX4RC.h"  // proprietary VHX4RC library for microscope control
#include <thread>
#include <chrono>
#include <math.h>
#include <vector>
#include <sstream>
#include <ctime>
#include "Eigen/Dense" 
#include <tuple>
#include <iomanip>
#include <cmath>

using namespace Eigen;

// same global varibles
bool bLongStep  = 1; // move in high speed mode
bool bSmallStep = 0; // move in slow speed mode    
int nLight = 6; // 6 = full coaxial
int nBrightnessValue = 120; // 120 == 200 in Auto scale
int nSize = 3;
int nImageModeTIFF = 1; // 1 == TIF
int nImageModeJPEG = 0; // 0 == JPEG

// lense powers
int nLensPower_x50 = 50;
int nLensPower_x700 = 700;
int nLensPower_x1500 = 1500;

// stability treshold
double RelaxationThreshold_x700 = 4.0; // threshold for stability check [um]
double RelaxationThreshold_x50 = 150.0; // threshold for stability check [um]

//------------------------------------------------------------------------------------------------------------------------------
// manual focus pos x50
double X_x50 = 13415; 
double Y_x50 = 8420;
double ZFocus_x50 = 18708.6; // lens position for the focus at x50

// grid postions 
double y_positions_x50[] = {-20080.0, -15060.0, -10040.0,  -5020.0, 0.0, 5020.0, 10040.0, 15060.0, 20080.0};
double x_positions_x50[] = {-19560.0, -13040.0,  -6520.0, 0.0, 6520.0,  13040.0,  19560.0};
double W_photo_x50 = 6.02*1000; // [um], width of photo (x50)
double H_photo_x50 = 4.52*1000; // [um], height of photo (x50)


// print in output the value of the variable VHX4RC_RESULT 
void printResult(HRESULT result) {
    switch (result)
        {
            case static_cast<HRESULT>(VHX4RC_S_OK):
                printf("VHX4RC_S_OK\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_COMMAND_FAILED):
                printf("VHX4RC_E_COMMAND_FAILED\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_ALREADY_INITIALIZED):
                printf("VHX4RC_E_ALREADY_INITIALIZED\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_NOT_INITIALIZED):
                printf("VHX4RC_E_NOT_INITIALIZED\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_COMM_ERROR):
                printf("VHX4RC_E_COMM_ERROR\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_INVALID_ARGUMENT):
                printf("VHX4RC_E_INVALID_ARGUMENT\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_OUTOFMEMORY):
                printf("VHX4RC_E_OUTOFMEMORY\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_PATH_NOT_FOUND):
                printf("VHX4RC_E_PATH_NOT_FOUND\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_SAVE_FILE_FAILED):
                printf("VHX4RC_E_SAVE_FILE_FAILED\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_FILE_ALREADY_EXISTS):
                printf("VHX4RC_E_FILE_ALREADY_EXISTS\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_INVALID_MODE):
                printf("VHX4RC_E_INVALID_MODE\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_INVALID_STATE):
                printf("VHX4RC_E_INVALID_STATE\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_Z_STAGE_NOTHING):
                printf("VHX4RC_E_Z_STAGE_NOTHING\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_Z_STAGE_DISCONNECTED):
                printf("VHX4RC_E_Z_STAGE_DISCONNECTED\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_INVALID_HANDLE):
                printf("VHX4RC_E_INVALID_HANDLE\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_NOT_SUPPORT):
                printf("VHX4RC_E_NOT_SUPPORT\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_OUTOFCOMPLIMIT):
                printf("VHX4RC_E_OUTOFCOMPLIMIT\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_XY_STAGE_DISCONNECTED):
                printf("VHX4RC_E_XY_STAGE_DISCONNECTED\n");
                break;
            case static_cast<HRESULT>(VHX4RC_E_XY_STAGE_NOT_INITIALIZED):
                printf("VHX4RC_E_XY_STAGE_NOT_INITIALIZED\n");
                break;
        }
}

// function to connect to the microscope 
void initConnection(){
    std::cout << "Program started..." << std::endl;
    // ------------- Connnect to microscope -------------
    const char* ipAddress = "XXX.XXX.XX.XXX"; // replace with the actual IP address of the VHX microscope
    std::cout << "\nConnecting to VHX microscope with IP=" << ipAddress << std::endl;
    HRESULT RESULT_InitA = VHX4RC_InitA(ipAddress);
    printResult(RESULT_InitA);
}

// function to initialize the XYZ plate
void InitializeXYZplate(){    
    std::cout  << "\nInitilize XYZ stage..." << std::endl;
    // declare two example variables of type HWND and DWORD
    HWND hWnd = GetConsoleWindow();   // Get handle to the console window for demonstration
    DWORD dwMessageId = WM_USER + 1;  // Use WM_USER + some value as an example
    HRESULT Result_StageInit = VHX4RC_XYStageInitializeThetaOrientation(hWnd, dwMessageId);
    printResult(Result_StageInit);
}

// function that moves the head to target position FinalPos.
// return the difference from the EndPosition to the target position
double MoveHeadTo(double FinalPos){
    double InitPos;
    double* pInitPos = &InitPos;

    VHX4RC_GetLensZPos(pInitPos); // get the initial lens position
    std::cout << "Starting: " << InitPos << std::endl;
    std::cout << "Target: " << FinalPos << std::endl;
    double deltaPos = 0.0;
    deltaPos = FinalPos - InitPos;
    if (deltaPos > 0) {
        while (deltaPos > 0) {
            if (deltaPos > 1) {
                VHX4RC_LensZUp(bLongStep);
            }
            else {
                VHX4RC_LensZUp(bSmallStep);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            VHX4RC_GetLensZPos(pInitPos);
            deltaPos = FinalPos - InitPos;
        }
    }
    else if (deltaPos < 0) {
        while (deltaPos < 0) {
            if (deltaPos < -1) {
                VHX4RC_LensZDown(bLongStep);
            }
            else {
                VHX4RC_LensZDown(bSmallStep);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            VHX4RC_GetLensZPos(pInitPos);
            deltaPos = FinalPos - InitPos;
        }
    }

    return deltaPos;
}

// set the illumation settings
void SetIllumination(int nLight, int nBrightnessValue){
    // set light to full
    std::cout << "\nSetting illumination..." << std::endl; 
    HRESULT Result_SetLight = VHX4RC_SetLightShift(nLight); 
    printResult(Result_SetLight);    
    HRESULT Result_SetBrightnessValue = VHX4RC_SetBrightnessValue(nBrightnessValue); 
    std::this_thread::sleep_for(std::chrono::seconds(5));
    printResult(Result_SetBrightnessValue);  
}


// ------------------------------------ FIT ------------------------------------------
// Function to fit a 2D quadratic function to the given points
VectorXd fit_quadratic(VectorXd x, VectorXd y, VectorXd z) {
    // Perform curve fitting
    VectorXd xy = (x.array() * y.array()).matrix();
    MatrixXd A(x.size(), 6);
    A << VectorXd::Ones(x.size()), x, y, xy, x.array().square(), y.array().square();
    VectorXd popt = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(z);
    return popt;
}

// function that creates the folder if not existing
void create_folder(const std::string folder){
    // check if the folder exists 
    bool isFolderExisting = std::filesystem::exists(folder);
    // create if it doesnt exist
    if (!isFolderExisting){
        std::filesystem::create_directories(folder);
    }
}


// function that performs the autofocus and returns the position of the lens
double AutoFocusProcedure(double nLastPosZ, int nLensPower, HWND hWnd,  DWORD dwMessageId) {    // Set lens power
    VHX4RC_SetLensPower(nLensPower);
    std::this_thread::sleep_for(std::chrono::seconds(3));

    // Start autofocus mode
    VHX4RC_StartAutoFocusMode();
    VHX4RC_SetAutoFocusAreaSize(nSize);
    std::this_thread::sleep_for(std::chrono::seconds(3));

    // Start the autofocus procedure
    HRESULT resultautofocus = VHX4RC_StartAutoFocus(hWnd, dwMessageId);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    printResult(resultautofocus);

    // Get the focus position
    double LensFinalPos;
    double* pLensFinalPos = &LensFinalPos;
    VHX4RC_GetLensZPos(pLensFinalPos);

    // If autofocus fails (difference > 50um) repeat until the difference is small or max attempts reached
    double focusDifference = LensFinalPos - nLastPosZ;
    double deltastep = -5.0;
    int attemptCount = 0;
    const int maxAttempts = 10;
    while (std::abs(focusDifference) > 50 && attemptCount <= maxAttempts) {
        double ztmp = MoveHeadTo(nLastPosZ + deltastep);
        VHX4RC_StartAutoFocus(hWnd, dwMessageId); // redo autofocus
        std::this_thread::sleep_for(std::chrono::seconds(20));
        VHX4RC_GetLensZPos(pLensFinalPos); // get head position
        // update difference
        focusDifference = LensFinalPos - nLastPosZ;
        deltastep = deltastep + 1.0;
        attemptCount++;
    }

    // Stop the autofocus mode
    VHX4RC_StopAutoFocusMode();

    // Return the final lens position
    return LensFinalPos;
}

// function that returns the formatted date and time
// main
int main() {
    std::cout << "Main started..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // photo save path
    std::string sample = "sample_name";
    std::string pathSavePhotos_check = "your_path_to_save_check"; // photos during the stabilization
    std::string pathSavePhotos = "your_path_to_save"; // real photos acquisition
    std::string pathFocusParams = "your_path_to_save_params"; // data for stability and fit

    // create the folders
    create_folder(pathSavePhotos_check);
    create_folder(pathSavePhotos);
    create_folder(pathFocusParams);    

    // connect to microscope
    initConnection();
    InitializeXYZplate();
    std::this_thread::sleep_for(std::chrono::seconds(5)); // wait for initialization to finish (time may be reduced)

    // set light to full
    SetIllumination(nLight, nBrightnessValue);
    std::this_thread::sleep_for(std::chrono::seconds(5));
    
    // get console hanlder and message id
    HWND hWnd = GetConsoleWindow();     // Get handle to the console window for demonstration 
    DWORD dwMessageId = WM_USER + 1;    // Use WM_USER + some value as an example

    // Data points used for focusing 
    int num_of_particles = 7; // number of particles/points
    VectorXd PosX(num_of_particles), PosY(num_of_particles), LastPosZ_x1500(num_of_particles), LastPosZ_x700(num_of_particles);
    // --- modify these vectors with your actual data ---
    PosX << 46, -26222, -33691, -33626, -33166, -23495, 34443;
    PosY <<  37328, 22727, 10435, 6428, 3096, -28652, -12193;
    LastPosZ_x1500 << 18645.8, 18640.6, 18643.1, 18645.2, 18647.9, 18645.3, 18663.7; //zPos manually registered
    LastPosZ_x700 = LastPosZ_x1500;

    // Define a vector of VectorXd to store multiple LastPosZ vectors --> memory of all LastPosZ to study stability
    std::vector<VectorXd> StorageZ_x1500;
    std::vector<VectorXd> StorageZ_x700;
    StorageZ_x1500.push_back(LastPosZ_x1500);
    StorageZ_x700.push_back(LastPosZ_x700);
    
    //----------------------- FIRST FIT WITH MANUAL MEASUREMENTS -----------------------------
    // Open the file to save parameters
    std::string fit_filename = pathFocusParams + "fit_parameters.txt";
    std::ofstream outfile(fit_filename, std::ios_base::app);

    // Write column names to the file
    outfile << "Date\tTime\ta\tb\tc\td\te\tf\n";

    // Print current time
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm* localTime = std::localtime(&now); // Convert the current time to local time
    
    // Format the date as "day/month"
    std::ostringstream dateStream;
    dateStream << std::setfill('0') << std::setw(2) << localTime->tm_mday << "/"; // Day
    dateStream << std::setfill('0') << std::setw(2) << (localTime->tm_mon + 1); // Month
    std::string formattedDate = dateStream.str();
    // Format the time as "hh:mm"
    std::ostringstream timeStream;
    timeStream << std::setfill('0') << std::setw(2) << localTime->tm_hour << ":"; // Hours
    timeStream << std::setfill('0') << std::setw(2) << localTime->tm_min; // Minutes
    std::string formattedTime = timeStream.str();

    //std::tie(formattedDate, formattedTime) = getFormattedDateTime(now);

    // First fit with a 2D quadratic function (fit performed on x700)
    VectorXd popt = fit_quadratic(PosX, PosY, LastPosZ_x700);

    // Append parameters to the file
    outfile << formattedDate << "\t" << formattedTime << "\t";
    for (int i = 0; i < 6; ++i) {
        outfile << popt(i) << "\t";
    }
    outfile << "\n";
    outfile.close();

    //----------------------- START WITH PROCEDURE (untill it stabilizes) -----------------------------

    // def fit parameters
    double a = 0;
    double b = 0;
    double c = 0;
    double d = 0;
    double e = 0;
    double f = 0;   
    double a_stable = 0;
    double b_stable = 0;
    double c_stable = 0;
    double d_stable = 0;
    double e_stable = 0;
    double f_stable = 0;
    auto startTime1 = std::chrono::steady_clock::now();

    int dataInterval = 15; //15 min
    double h_max = 6; //max num of hours for the code to trying to reach stability
    int numOfiterations_x700 = (h_max * 60)/(dataInterval); //max num of iterations for the x700 loop

    // stabiliy check
    int counter_x700 = 0; //counter for the x700 loop
    bool isStable_x700 = false; // Flag to track stability for the x700 loop
    while (!isStable_x700 && counter_x700 <= numOfiterations_x700) {
        // Get current time
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::minutes>(currentTime - startTime1);
        
        // If 15 minutes have elapsed, print current time and reset start time
        if (elapsedTime.count() >= dataInterval || counter_x700 == 0) {
            // Print current time
            auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << "Current time: " << std::ctime(&now);
            
            // loop over particles
            for (int i = 0; i < PosX.size(); i++) {

                // 1) get the xyz positions 
                int nPosX = PosX[i];
                int nPosY = PosY[i];
                double nLastPosZ_x700 = LastPosZ_x700[i];                
                double nLastPosZ_x1500 = LastPosZ_x1500[i];                
                std::string stability_filename = pathFocusParams + "stability_check_time_focus_part_" + std::to_string(i) + ".txt";

                // 2) move to x, y position
                VHX4RC_XYStageMoveTo(nPosX, nPosY);
                std::cout << "nPosX = " << nPosX << std::endl;
                std::cout << "nPosY = " << nPosY << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(10)); // wait for x,y posituion to be reached

                // 3) move to last i-th focus height 
                double z = MoveHeadTo(nLastPosZ_x1500);  
                std::this_thread::sleep_for(std::chrono::seconds(5));

                // 4) ------ autofocus at x1500 ------
                double LensFocusPos_x1500 = AutoFocusProcedure(nLastPosZ_x700, nLensPower_x1500, hWnd, dwMessageId);

                // 5) ------ autofocus at x700 ------
                double LensFocusPos_x700 = AutoFocusProcedure(LensFocusPos_x1500, nLensPower_x700, hWnd, dwMessageId);

                // put last focus position in i-th positino of z array
                LastPosZ_x700[i] = LensFocusPos_x700;
                LastPosZ_x1500[i] = LensFocusPos_x1500;

                // Format the date as "day/month"
                
                std::tm* localTime = std::localtime(&now); // Convert the current time to local time
                std::ostringstream dateStream;
                dateStream << std::setfill('0') << std::setw(2) << localTime->tm_mday << "/"; // Day
                dateStream << std::setfill('0') << std::setw(2) << (localTime->tm_mon + 1); // Month
                std::string formattedDate = dateStream.str();
                // Format the time as "hh:mm"
                std::ostringstream timeStream;
                timeStream << std::setfill('0') << std::setw(2) << localTime->tm_hour << ":"; // Hours
                timeStream << std::setfill('0') << std::setw(2) << localTime->tm_min; // Minutes
                std::string formattedTime = timeStream.str();
                //std::tie(formattedDate, formattedTime) = getFormattedDateTime(now);

                // Save time and focus position in a txt file
                std::ofstream file(stability_filename, std::ios_base::app);
                if (file.is_open()) {
                    // Write the formatted time to the file and the focus position
                    file << formattedDate << "\t" << formattedTime << "\t" << LensFocusPos_x1500 << "\t" << LensFocusPos_x700 << std::endl;
                    file.close();
                } else {
                    std::cerr << "Error: Unable to open file!" << std::endl;
                };

                // take photo (photos are to check if the focus is correct, not used for the fit)
                std::string file_stability_check = "time_focus_part_" + std::to_string(i);
                std::string fileName = pathSavePhotos_check + file_stability_check + "_" + std::to_string(counter_x700) + "_z" + std::to_string(LensFocusPos_x700) + ".tif";
                const char* filenamePtr = fileName.c_str();
                HRESULT Result_Saveimage = VHX4RC_SaveImageVHXHDDA(nImageModeTIFF, filenamePtr);
                std::this_thread::sleep_for(std::chrono::seconds(5));
            }
            StorageZ_x700.push_back(LastPosZ_x700);
            StorageZ_x1500.push_back(LastPosZ_x1500);

            // --------------------- FIT -------------------

            // Convert the current time to local time
            std::tm* localTime = std::localtime(&now);
            // Format the date as "day/month"
            std::ostringstream dateStream;
            dateStream << std::setfill('0') << std::setw(2) << localTime->tm_mday << "/"; // Day
            dateStream << std::setfill('0') << std::setw(2) << (localTime->tm_mon + 1); // Month
            std::string formattedDate = dateStream.str();
            // Format the time as "hh:mm"
            std::ostringstream timeStream;
            timeStream << std::setfill('0') << std::setw(2) << localTime->tm_hour << ":"; // Hours
            timeStream << std::setfill('0') << std::setw(2) << localTime->tm_min; // Minutes
            std::string formattedTime = timeStream.str(); 
            //std::tie(formattedDate, formattedTime) = getFormattedDateTime(now);

            // Fit the height points with a 2D quadratic function
            VectorXd popt = fit_quadratic(PosX, PosY, LastPosZ_x700);
            
            // save parameters to the file  
            std::ofstream outfile(fit_filename, std::ios_base::app);
            outfile << formattedDate << "\t" << formattedTime << "\t";
            for (int i = 0; i < 6; ++i) {
                outfile << popt(i) << "\t";
            }
            outfile << "\n";
            outfile.close();

            // update fit parameters
            a = popt(0);
            b = popt(1);
            c = popt(2);
            d = popt(3);
            e = popt(4);
            f = popt(5);

            // Reset start time and update the loop
            startTime1 = currentTime;            
            counter_x700 = counter_x700 + 1;

            // Check for stability
            if (StorageZ_x700.size() < 5) {
                isStable_x700 = false;

            } else {
                // Vectors for the last three zPos
                int start_index = StorageZ_x700.size() - 5; 
                VectorXd z_1 = StorageZ_x700[start_index];
                VectorXd z_2 = StorageZ_x700[start_index + 1];
                VectorXd z_3 = StorageZ_x700[start_index + 2];
                VectorXd z_4 = StorageZ_x700[start_index + 3];
                VectorXd z_5 = StorageZ_x700[start_index + 4];
                
                // Flag to track if all differences are less than 1 um
                bool allDifferencesLessThanOne = true; 


                // Updated stablity check: within 5 points, z_max - z_min < 1um for each particle
                for (int ii = 0; ii < num_of_particles; ii++) {

                    VectorXd z_pos_ith(5); // vector to store last 5 z_pos for each particles
                    double z_max = 0;
                    double z_min = 200000;
                    double z_disp = 0;
                    
                    z_pos_ith[0] = z_1[ii];
                    z_pos_ith[1] = z_2[ii];
                    z_pos_ith[2] = z_3[ii];
                    z_pos_ith[3] = z_4[ii];
                    z_pos_ith[4] = z_5[ii];

                    for (int j = 0; j < z_pos_ith.size(); j++) {
                        if (z_pos_ith[j] > z_max) {
                            z_max = z_pos_ith[j];
                        }

                        if (z_pos_ith[j] < z_min) {
                            z_min = z_pos_ith[j];
                        }
                    }

                    z_disp = z_max - z_min;
                    std::cout << "Particle " << ii << ": z_max = " << z_max << ", z_min = " << z_min << ", z_disp = " << z_disp << std::endl;

                    if (z_disp >= RelaxationThreshold_x700) {
                        allDifferencesLessThanOne = false; // Set flag to false if any z_max - z_min is not less than 1 um
                        break; // Exit the loop as soon as one difference is not less than 1 um
                    }                
                }

                // If all differences are less than 1 um for every particle, declare stability
                if (allDifferencesLessThanOne) {
                    isStable_x700 = true;
                    std::cout << "Stability x700 reached." << std::endl;

                    // fit performed on the average of the last 5 z_positions for each particle
                    VectorXd z_average(num_of_particles);
                    for (int jj = 0; jj < num_of_particles; jj++) {
                        z_average[jj] = (z_1[jj] + z_2[jj] + z_3[jj] + z_4[jj] + z_5[jj]) / 5;
                    }

                    // Repeat the correct fit
                    VectorXd popt_stable = fit_quadratic(PosX, PosY, z_average);
                    a_stable = popt_stable(0);
                    b_stable = popt_stable(1);
                    c_stable = popt_stable(2);
                    d_stable = popt_stable(3);
                    e_stable = popt_stable(4);
                    f_stable = popt_stable(5);
                }
            }
        }
    }

    /*
    at a certain point you either exit because:
    - stability reached
    - max number of iterations reached
    */
    if (counter_x700 > numOfiterations_x700) {

        // stop the autofocus mode
        HRESULT resultstopautofocusmode = VHX4RC_StopAutoFocusMode();
        printResult(resultstopautofocusmode);

        // disconnect from the microscope and stop
        std::cout << "\nDisconnection..." << std::endl; 
        bool disconnection_bool = 0;
        HRESULT Result_exit = VHX4RC_Exit(disconnection_bool);
        printResult(Result_exit); 

        // Output error message
        std::cout << "Stability x700 not reached." << std::endl;

        return 0;   
    }
    
    // Sleep for a short duration before checking again
    std::this_thread::sleep_for(std::chrono::seconds(1));

    // --------------------------------- CONTINUE WIHT PROCEDURE (photos) ---------------------------------
    // -------------------------------------- x700 -------------------------------------------------------

    double LensFinalPos_x700;
    double* pLensFinalPos_x700 = &LensFinalPos_x700;

    for (const auto &x : x_positions_x50) {
        for (const auto &y : y_positions_x50) {
            double deltaX = W_photo_x50/5;
            double deltaY = H_photo_x50/4;
            double x_positions_x700[] = {x-2*deltaX, x-1*deltaX, x, x+1*deltaX, x+2*deltaX};  
            double y_positions_x700[] = {y-1*deltaY, y, y+1*deltaY};
            // y --> H
            // x --> W


            // loop over the sub grid of photos for x700
            for (const auto &nPosX : x_positions_x700) {
                for (const auto &nPosY : y_positions_x700) {

                    if (nPosX == x && nPosY != y) {
                        // skip this position
                        continue;
                    }
                    else {                                                
                        // move head to focus position
                        double z = a_stable + b_stable*nPosX + c_stable*nPosY + d_stable*nPosX*nPosY + e_stable*nPosX*nPosX + f_stable*nPosY*nPosY; // calculate the z position
                        double res = MoveHeadTo(z); // move head   
                        std::this_thread::sleep_for(std::chrono::milliseconds(500));

                        // move the plate to XY position
                        HRESULT Result_XYStageMoveTo = VHX4RC_XYStageMoveTo(nPosX, nPosY); 
                        std::cout << "   x=" << nPosX << "; y=" << nPosY << std::endl;
                        
                        // wait 5s only if moving to a new grid (must move for longer) else wait 0.5s (moving inside small grid takes less)
                        bool IsFirstPhotoOfNewSubgrid = ((nPosX == x_positions_x700[0]) && (nPosY == y_positions_x700[0]));
                        if ( IsFirstPhotoOfNewSubgrid == true ) {
                            std::this_thread::sleep_for(std::chrono::seconds(5));
                            }
                        else {
                        std::this_thread::sleep_for(std::chrono::milliseconds(500));
                            }

                        // get lens position
                        VHX4RC_GetLensZPos(pLensFinalPos_x700);

                        // take photo (set on .tif)
                        std::string fileName = pathSavePhotos + sample + "_X700" + "_x" + std::to_string(nPosX) + "_y" + std::to_string(nPosY) + "_z" + std::to_string(LensFinalPos_x700) + ".tif";
                        const char* filenamePtr = fileName.c_str();
                        HRESULT Result_Saveimage = VHX4RC_SaveImageVHXHDDA(nImageModeTIFF, filenamePtr);
                        printResult(Result_Saveimage);
                    }
                }
            }
        }
    }
    std::this_thread::sleep_for(std::chrono::seconds(5));

    // -------------------------------------- x50 ---------------------------------------------------------

    // set lens power
    VHX4RC_SetLensPower(nLensPower_x50);
    std::cout << "lens: " << nLensPower_x50 << std::endl; 
    std::this_thread::sleep_for(std::chrono::seconds(5));

    // Move to x, y position
    VHX4RC_XYStageMoveTo(X_x50, Y_x50);

    // Move head to the manually measured value: ZFocus_x50 
    double restore = MoveHeadTo(ZFocus_x50); // move head   
    std::this_thread::sleep_for(std::chrono::seconds(10));

    bool isStable_x50 = false; // Flag to track stability
    int counter_x50 = 0;
    int numOfiterations_x50 = 10; //max num of iterations for the x50 loop --> fixed to 10

    /*
    while NOT-STABLE and STILL-IN-COUNTS do procedure of autofocusing and check for stability
    */
    while (!isStable_x50 && counter_x50 <= numOfiterations_x50) {
        // autofocus to check the focus at x50
        // start autofocus mode
        HRESULT resultstartautofocusmode = VHX4RC_StartAutoFocusMode();
        HRESULT resultsetregionsize = VHX4RC_SetAutoFocusAreaSize(nSize);
        std::this_thread::sleep_for(std::chrono::seconds(5));

        // start the autofocus
        HRESULT resultautofocus = VHX4RC_StartAutoFocus(hWnd, dwMessageId);
        std::this_thread::sleep_for(std::chrono::seconds(20)); // wait for focus determination to finish
        printResult(resultautofocus);              
        
        // Stop the autofocus mode
        HRESULT resultstopautofocusmode = VHX4RC_StopAutoFocusMode();
        printResult(resultstopautofocusmode);

        // get the focus position
        double LensFinalPos_x50;
        double* pLensFinalPos_x50 = &LensFinalPos_x50;
        VHX4RC_GetLensZPos(pLensFinalPos_x50);

        std::cout << "z_pos_x50 = " << LensFinalPos_x50 << std::endl;        

        double diff = ZFocus_x50 - LensFinalPos_x50; // difference between the autofocus value and the manually given position
        if (std::abs(diff) < RelaxationThreshold_x50) {
            // set stability bool to true
            isStable_x50 = true;
            std::cout << "Stability x50 reached." << std::endl;
            double restore = MoveHeadTo(LensFinalPos_x50); // move head   

            // Take photos if the autofocus value is within 150 of the manually given position
            for (const auto &nPosX : x_positions_x50) {
                for (const auto &nPosY : y_positions_x50) {

                    // move the plate to XY position
                    HRESULT Result_XYStageMoveTo = VHX4RC_XYStageMoveTo(nPosX, nPosY); 
                    std::cout << "x=" << nPosX << "; y=" << nPosY << std::endl;
                    std::this_thread::sleep_for(std::chrono::seconds(5));

                    // take photo (set on .jpg)
                    std::string fileName = pathSavePhotos + sample + "_X50" + "_x" + std::to_string(nPosX) + "_y" + std::to_string(nPosY) + "_z" + std::to_string(LensFinalPos_x50) + ".tif";                    
                    const char* filenamePtr = fileName.c_str();
                    HRESULT Result_Saveimage = VHX4RC_SaveImageVHXHDDA(nImageModeTIFF, filenamePtr);
                    printResult(Result_Saveimage);
                    }
                }

            } else {
            // If autofocus doesn't converge or differs by more than 100, increment the counter
            counter_x50++;
            }
        
    }
    if (counter_x50 > numOfiterations_x50) {
        // Output error message
        std::cout << "Stability x50 not reached." << std::endl;
        }


    // ----------------------------------- DISCONNTECT ------------------------------------------------------
    // disconnect from microscope
    std::cout << "\nDisconnection..." << std::endl; 
    bool disconnection_bool = 0;
    HRESULT Result_exit = VHX4RC_Exit(disconnection_bool);
    printResult(Result_exit); 

    // stop the timer
    auto end      = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    int hour      = std::chrono::duration_cast<std::chrono::hours>(duration).count();
    duration     -= std::chrono::hours(hour);
    int minutes   = std::chrono::duration_cast<std::chrono::minutes>(duration).count();
    duration     -= std::chrono::minutes(minutes);
    int seconds   = duration.count();
 
    // Output the duration
    std::cout << "Time taken by the process: " << hour << " hours, " << minutes << " minutes, and " << seconds << " seconds" << std::endl;
    return 0;   

}