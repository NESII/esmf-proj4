/* 
 * File:   test_transform.cpp
 * Author: ben.koziol
 *
 * Created on Apr 7, 2015, 11:16:06 AM
 */

#include <iostream>
#include <sstream>
#include "transform.h"
#include "test_transform_data.h"

using namespace std;

// #####################################################################################################################
// Test-suite globals.

// precision of standard out used for printing
const int COUT_PRECISION = 16;
// TODO is this tolerance setting okay?
// tolerance setting for checking value assertions. if the difference is less than the tolerance, tests will pass.
const double TEST_TOLERANCE = 1e-12;
int TEST_TOTAL;
int TEST_SKIPPED;
int TEST_FAILED;

// #####################################################################################################################
// Test-related PROJ4 strings.

string PROJ_SPHERICAL = "+proj=longlat +a=6370997 +b=6370997 +towgs84=0,0,0,0,0,0,0 +no_defs ";
string PROJ_LCC = "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=47.5 +lon_0=-97 +x_0=3325000 +y_0=2700000 +ellps=WGS84 +units=m +no_defs ";
string PROJ_WGS84 = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs ";
string PROJ_SPHERICAL_OVER = "+proj=longlat +a=6370997 +b=6370997 +towgs84=0,0,0,0,0,0,0 +no_defs +lon_wrap=180 ";
string PROJ_WGS84_OVER = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs +lon_wrap=180 ";
string PROJ_SPHERICAL_SHORT = "+proj=longlat +a=6370997 +b=6370997 ";
string PROJ_WGS84_SHORT = "+proj=longlat +ellps=WGS84 ";

// #####################################################################################################################

class AssertionError: public exception{
    int line;
    public:
        AssertionError(){this->line = -1;};
        AssertionError(int line){this->line = line;};
        virtual const char* what() const throw(){
            stringstream str;
            str << "AssertionError: line " << this->line << endl;
            const char * result = str.str().c_str();
            return result;
        }
};


class SkipTest: public exception{
    const char* msg;
    public:
        SkipTest(const char* msg){this->msg = msg;};
        virtual const char* what() const throw(){
            stringstream str;
            str << "SkipTest: " << this->msg << endl;
            const char* result = str.str().c_str();
            return result;
        }
};


void assertDoubleEqual(double a, double b, int line = -1){
    if (a != b){
        AssertionError bad_assertion(line);
        throw bad_assertion;
    }
}


void assertEqual(double a, double b, int line = -1){
    if (a != b){
        AssertionError bad_assertion(line);
        throw bad_assertion;
    }
}

void assertGeneral(const bool& result, const int& line = -1){
    if (!result){
        AssertionError ae(line);
        throw ae;
    }
}

void assertTrue(const bool& result, const int& line = -1){
    if (!result){
        AssertionError ae(line);
        throw ae;
    }
}

void assertFalse(const bool& result, const int& line = -1){
    if (result){
        AssertionError ae(line);
        throw ae;
    }
}

void assertVectorEqual(double vec[], double vec_actual[], int size, int line = -1, double tolerance = TEST_TOLERANCE){
    for (int i = 0; i < size; i++) {
        if (vec[i] != vec_actual[i]){
            double difference = vec[i] - vec_actual[i];
            if (abs(difference) > tolerance){
                cout.precision(COUT_PRECISION);
                cout << "index: " << i << endl;
                cout << "vec: " << vec[i] << endl;
                cout << "vec_actual: " << vec_actual[i] << endl;
                cout << "difference (vec - vec_actual): " << difference << endl;
                AssertionError bad_assertion(line);
                throw bad_assertion;
            }
        }
    }
}

void fillVector(double src[], int size, double dst[]){
    for (int i = 0; i < size; i++){
        dst[i] = src[i];
    }
}

void interleave2d(double vec[], int size, double vec_a[], double vec_b[]){
    for (int i = 0; i < size; i += 2){
        vec[i] = vec_a[i/2];
        vec[i+1] = vec_b[i/2];
    }
}

void printVector(double vec[], int size){
    for (int i = 0; i < size; i++){
        cout << vec[i] << endl;
    }
}

void vectorMultiply(double vec[], int size, double b){
    for (int i = 0; i < size; i++){
        vec[i] *= b;
    }
}

// #####################################################################################################################
// #####################################################################################################################

class TestData: private Uncopyable {
    public:
        Coordinates* p_coords;
        TestData();
        TestData(double values_degrees[], int size, int x_index, int y_index);
        TestData(double values_degrees[], int size, int x_index, int y_index, int z_index);
		TestData(double values_degrees[], int size, int x_index, int y_index, int z_index, CoordinateSystem& coordsys);
        ~TestData();
	private:
		CoordinateSystem* p_coordsys;
		CoordinateIndex* p_coord_index;
        double* p_values;
        bool has_new_values;
        void initialize(double values_degrees[], int size, int x_index, int y_index, int z_index, bool has_new_values, CoordinateSystem* coordsys);
};

TestData::TestData() {
    int size = 6;
    double* values = new double[size];
    values[0] = 75;
    values[1] = -170;
    values[2] = 80;
    values[3] = -160;
    values[4] = 85;
    values[5] = -150;
    this->p_values = values;
    this->initialize(values, size, 1, 0, -1, true, NULL);
};

TestData::TestData(double values_degrees[], int size, int x_index, int y_index){
    this->initialize(values_degrees, size, x_index, y_index, -1, false, NULL);
};

TestData::TestData(double values_degrees[], int size, int x_index, int y_index, int z_index){
    this->initialize(values_degrees, size, x_index, y_index, z_index, false, NULL);
};

TestData::TestData(double values_degrees[], int size, int x_index, int y_index, int z_index, CoordinateSystem& coordsys){
    this->initialize(values_degrees, size, x_index, y_index, z_index, false, &coordsys);
};

TestData::~TestData() {
    if (this->has_new_values){
        delete[] this->p_values;
    }
	delete this->p_coord_index;
	if (this->p_coordsys != NULL){
		delete this->p_coordsys;
	}
    delete this->p_coords;
};

void TestData::initialize(double values_degrees[], int size, int x_index, int y_index, int z_index, bool has_new_values,
                          CoordinateSystem* coordsys){
    this->has_new_values = has_new_values;
	vector<int> z_exclude;
	bool skip = false;

	if (z_index != -1){
		for (int i = z_index; i < size; i+=3){
			z_exclude.push_back(i);
		}
	}

    for (int i = 0; i < size; i++){
		skip = false;
		if (z_index != -1){
			for (vector<int>::const_iterator it = z_exclude.begin(); it != z_exclude.end(); ++it){
				if (*it == i){
					skip = true;
					break;
				}
			}
		}
		if (skip == false){
			values_degrees[i] *= DEG_TO_RAD;
		}
    }

	this->p_coord_index = new CoordinateIndex(x_index, y_index, z_index);
	
	if (coordsys == NULL){
		this->p_coordsys = new CoordinateSystem(PROJ_SPHERICAL);
		coordsys = this->p_coordsys;
	}
    else {
		this->p_coordsys = NULL;
	}
	
    Coordinates* coords = new Coordinates(values_degrees, size, *this->p_coord_index, *coordsys);
    this->p_coords = coords;
};

// #####################################################################################################################
// #####################################################################################################################

void test_Coordinates(){
    TestData td;
    Coordinates *p_coords = td.p_coords;
    const CoordinateIndex* cindex = p_coords->getCoordinateIndexPtr();
    double *p_values = p_coords->getValuesPtr();
    assertTrue(p_values[cindex->getXLocation()]*RAD_TO_DEG == -170, __LINE__);
    assertTrue(p_values[cindex->getYLocation()]*RAD_TO_DEG == 75, __LINE__);
    assertTrue(p_coords->getNdim() == 2, __LINE__);
    assertTrue(p_coords->getSize() == 6, __LINE__);
    assertTrue(p_coords->getCoordinateSystemPtr()->getWrappedState() == WRAPPED, __LINE__);
};

void test_Coordinates_updateWrappedState(){
    TestData td;
    assertTrue(td.p_coords->getCoordinateSystemPtr()->getWrappedState() == WRAPPED, __LINE__);
};

void test_Coordinates_transformTo(){
    TestData td;
    double values_actual[] = {75.0957, -170, 80.0655, -160, 85.0332, -150};
    string proj_string_dest = "+proj=longlat +ellps=WGS84";
    CoordinateSystem coordsys_dest(proj_string_dest);
    td.p_coords->transformTo(coordsys_dest);
    assertTrue(*td.p_coords->getCoordinateSystemPtr() == coordsys_dest, __LINE__);

    double *p_values = td.p_coords->getValuesPtr();
    for (int i = 0; i < td.p_coords->getSize(); i++){
        p_values[i] *= RAD_TO_DEG;
    }

    assertVectorEqual(p_values, values_actual, td.p_coords->getSize(), __LINE__);
};

void test_Coordinates_transformTo_rotated_pole(){
    int size = 4;
    double values[] = {-25.075000762939453, -20.075000762939453, -24.96500015258789, -19.96500015258789};
    double values_actual[] = {-8.33782, 26.2072, -8.27533, 26.3473};
    string proj_dest = "+proj=ob_tran +o_proj=latlon +o_lon_p=-162.0 +o_lat_p=39.25 +lon_0=180 +ellps=sphere";

    // Convert the coordinate values to radians.
    for (int i = 0; i < size; i++){
        values[i] *= DEG_TO_RAD;
    }
    CoordinateSystem cs_src(PROJ_SPHERICAL);
    CoordinateIndex ci(0, 1, -1);
    Coordinates coords(values, 4, ci, cs_src);

    CoordinateSystem cs_dest(proj_dest);

    coords.transformTo(cs_dest);
    double *new_values = coords.getValuesPtr();
    // Convert the coordinate values back to degrees.
    for (int i = 0; i < size; i++){
        new_values[i] *= RAD_TO_DEG;
    }
    assertVectorEqual(new_values, values_actual, size);
};

void test_Coordinates_transformTo_various_transformations() {
    int size = 8;
    int ndim = 2;
    double coords_original[] = {-105.4, 39.4, -105.2, 39.4, -105.4, 39.8, -105.2, 39.8};
    double coords[] = {-105.4, 39.4, -105.2, 39.4, -105.4, 39.8, -105.2, 39.8};

    double x_lcc_actual[] = {2623434.05872470699250698089599609375, 2640078.8396765361540019512176513671875, \
                                   2627950.5938705527223646640777587890625, 2644488.219195603393018245697021484375};
    double y_lcc_actual[] = {1866370.6100018967, 1864638.7734324972, 1909258.7752076953, 1907538.087840777};
    double coords_lcc_actual[size];
    interleave2d(coords_lcc_actual, size, x_lcc_actual, y_lcc_actual);

    // transform coordinates from spherical to lambert conformal conic
    TestData td(coords, size, 0, 1, -1);
    CoordinateSystem coordsys_dest(PROJ_LCC);
    td.p_coords->transformTo(coordsys_dest);
    double *p_values = td.p_coords->getValuesPtr();
    assertVectorEqual(p_values, coords_lcc_actual, size, __LINE__);

    // transform coordinates from lambert conformal conic to spherical
    CoordinateSystem proj_spherical(PROJ_SPHERICAL);
    int pj_code = td.p_coords->transformTo(proj_spherical);
    assertTrue(pj_code == 0, __LINE__);
    vectorMultiply(coords, size, RAD_TO_DEG);
    assertVectorEqual(p_values, coords_original, size, __LINE__);

    // transform coordinates from wgs84 to spherical
    TestData td_roundtrip(coords, size, 0, 1, -1);
    CoordinateSystem proj_wgs84(PROJ_WGS84);
    td_roundtrip.p_coords->transformTo(proj_wgs84);
    td_roundtrip.p_coords->transformTo(proj_spherical);
    double *p_values_roundtrip = td_roundtrip.p_coords->getValuesPtr();
    vectorMultiply(p_values_roundtrip, size, RAD_TO_DEG);
    assertVectorEqual(p_values_roundtrip, coords_original, size, __LINE__);
};

void test_getWrappedState(){
    int size = 6;

    double coords[] = {0, 20, 5, 30, 10, 40};
    TestData td(coords, size, 0, 1);
    WrappedState ws = getWrappedState(*td.p_coords);
    assertGeneral(ws == UNDEFINED, __LINE__);

    double coords2[] = {290, 0, 30, 5, 40, 10};
    TestData td2(coords2, size, 1, 0);
    ws = getWrappedState(*td2.p_coords);
    assertGeneral(ws == UNDEFINED, __LINE__);

    double coords3[] = {290, 0, 30, 5, 40, 10};
    TestData td3(coords3, size, 0, 1);
    ws = getWrappedState(*td3.p_coords);
    assertGeneral(ws == UNWRAPPED, __LINE__);

    double coords4[] = {30, 0, -170, 5, 40, 10};
    TestData td4(coords4, size, 0, 1);
    ws = getWrappedState(*td4.p_coords);
    assertGeneral(ws == WRAPPED, __LINE__);

    size = 9;
    double coords5[] = {30, 0, 50, -170, 5, 40, 10, 80, 6};
    TestData td5(coords5, size, 0, 1, 2);
    ws = getWrappedState(*td5.p_coords);
    assertGeneral(ws == WRAPPED, __LINE__);
}

void test_CoordinateSystem(){
    // Test the towgs84 datum conversion is appropriately added.
    CoordinateSystem ep(PROJ_SPHERICAL_SHORT);
    string proj_string_actual = "+proj=longlat +a=6370997 +b=6370997 +towgs84=0,0,0";
    CoordinateSystem ep_actual(proj_string_actual);
    assertTrue(ep.getWrappedState() == UNINITIALIZED, __LINE__);
    assertTrue(ep == ep_actual, __LINE__);
    assertTrue(ep.isSphericalLongLat(), __LINE__);
    assertTrue(ep.isWrappableLongLat(), __LINE__);

    // Test initializing with a wrapped state.
    CoordinateSystem ep_wrap(PROJ_SPHERICAL_SHORT, WRAPPED);
    assertTrue(ep_wrap.getWrappedState() == WRAPPED, __LINE__);

    // Test no_defs parameter is added by default.
    if (ep.getValuePtr()->find("no_defs") == ep.getValuePtr()->end()) {
        AssertionError ae(__LINE__);
        throw ae;
    }

    // Test initializing with an unwrapped state.
    CoordinateSystem ep_unwrap(PROJ_SPHERICAL_SHORT, UNWRAPPED);
    CoordinateSystem ep_unwrap_actual("+a=6370997 +b=6370997 +lon_wrap=180 +proj=longlat +towgs84=0,0,0");
    assertTrue(ep_unwrap == ep_unwrap_actual, __LINE__);
    assertTrue(ep_unwrap.getWrappedState() == UNWRAPPED, __LINE__);

    // Test the wgs84 datum conversion is not added if present.
    CoordinateSystem ep2(PROJ_WGS84);
    bool result2 = ep2.isSphericalLongLat();
    assertFalse(result2, __LINE__);
    bool result2a = ep2.isWGS84LongLat();
    assertTrue(result2a, __LINE__);
    // Test this is a longitude/latitude coordinate system.
    bool result3 = ep2.isLongLat();
    assertTrue(result3);

    // Test a cartesian coordinate system is not longitude/latitude.
    CoordinateSystem ep3(PROJ_LCC);
    assertTrue(ep3 != ep);
    assertFalse(ep3.containsToWGS84(), __LINE__);
    assertFalse(ep3.isSphericalLongLat(), __LINE__);
    assertFalse(ep3.isWGS84LongLat(), __LINE__);
    assertFalse(ep3.isLongLat(), __LINE__);
    assertFalse(ep3.isWrappableLongLat(), __LINE__);

    // Test passing through without validation and updating.
    string proj_string_unvalidated = "+proj=longlat +ellps=WGS84";
    CoordinateSystem ep_no_validation(proj_string_unvalidated, UNWRAPPED, false);
    assertTrue(proj4MapToString(*ep_no_validation.getValuePtr()) == "+ellps=WGS84 +proj=longlat ", __LINE__);

    // Test setting the wrapped state.
    CoordinateSystem set_wrapping(proj_string_unvalidated);
    set_wrapping.setWrappedState(UNWRAPPED);
    assertTrue(set_wrapping.getWrappedState() == UNWRAPPED, __LINE__);
    CoordinateSystem set_wrapping_actual("+ellps=WGS84 +lon_wrap=180 +no_defs +proj=longlat +towgs84=0,0,0");
    assertTrue(set_wrapping == set_wrapping_actual, __LINE__);
	
	// Test wrapped state appropriately read with "lon_wrap" key.
	CoordinateSystem wrapped_coordsys("+proj=longlat +ellps=WGS84 +lon_wrap=180");
	assertTrue(wrapped_coordsys.getWrappedState() == UNWRAPPED, __LINE__);
}

void test_Coordinates_transformTo_bigger_arrays(){
    int size = 128*64*2;
    double values[size];
    double values_original[size];
    interleave2d(values, size, test_x_big, test_y_big);
    // TODO add actual (correct) coordinates
    fillVector(values_original, size, values);
    vectorMultiply(values, size, DEG_TO_RAD);

    CoordinateSystem proj_wgs84(PROJ_WGS84);
    CoordinateSystem proj_wgs84_over(PROJ_WGS84_OVER);
    CoordinateSystem proj_spherical(PROJ_SPHERICAL);
    CoordinateSystem proj_spherical_over(PROJ_SPHERICAL_OVER);

	TestData td(values, size, 0, 1, -1, proj_wgs84);
    td.p_coords->transformTo(proj_spherical_over);
    td.p_coords->transformTo(proj_wgs84_over);
    double *cvalues = td.p_coords->getValuesPtr();
    vectorMultiply(cvalues, size, RAD_TO_DEG);

    for (int i = 0; i < size; i += 2){
        if (cvalues[i] < 0){
            AssertionError ae(__LINE__);
            throw ae;
        }
    }

    assertVectorEqual(cvalues, values_original, size, __LINE__);
}

void test_Coordinates_transformTo_edge_effects(){
    coord_overloads cmap;
    int size_xy = 2;
    double x[] = {-180, 180};
    double y[] = {-90, 90};
    int size = 2*2*2;
    double coords[size];
    double coords_original[size];
    int ctr = 0;
    for (int i = 0; i < size_xy; i++){
        for (int j = 0; j < size_xy; j++){
            coords[ctr] = x[i];
            coords[ctr+1] = y[i];
            ctr += 2;
        }
    }
    fillVector(coords, size, coords_original);
    CoordinateSystem proj_src(PROJ_WGS84);
    CoordinateSystem proj_dst(PROJ_SPHERICAL);
    TestData td(coords, size, 0, 1, -1, proj_src);
    td.p_coords->transformTo(proj_dst);
    double *p_values = td.p_coords->getValuesPtr();
    vectorMultiply(p_values, size, RAD_TO_DEG);
    assertVectorEqual(coords, coords_original, size, __LINE__);
}

void test_Coordinates_transformTo_edge_effects_wrapped(){
    coord_overloads cmap;
    int size_xy = 2;
    double x[] = {0, 360};
    double y[] = {-90, 90};
    int size = 2*2*2;
    double coords[size];
    double coords_original[size];
    int ctr = 0;
    for (int i = 0; i < size_xy; i++){
        for (int j = 0; j < size_xy; j++){
            coords[ctr] = x[i];
            coords[ctr+1] = y[i];
            ctr += 2;
        }
    }

    fillVector(coords, size, coords_original);
    CoordinateSystem proj_src(PROJ_WGS84);
    CoordinateSystem proj_dst(PROJ_SPHERICAL_OVER);
    TestData td(coords, size, 0, 1, -1, proj_src);
    td.p_coords->transformTo(proj_dst);
    double *p_values = td.p_coords->getValuesPtr();
    vectorMultiply(p_values, size, RAD_TO_DEG);
    assertVectorEqual(coords, coords_original, size, __LINE__);
}

void test_Coordinates_transformTo_three_dimensions(){
    int size = 3;
    double coords_values[] = {-180, 40, 50};
    double coords_actual[] = {-180, 40, 0.049};
    CoordinateSystem proj_src("+proj=longlat +a=6370997 +b=6370997 +towgs84=0,0,0,0,0,0,0 +no_defs +vunits=m");
    CoordinateSystem proj_dst("+proj=longlat +a=6370998 +b=6370998 +towgs84=0,0,0,0,0,0,0 +no_defs +vunits=km");
    TestData td(coords_values, size, 0, 1, 2, proj_src);
    td.p_coords->transformTo(proj_dst);

    double *values = td.p_coords->getValuesPtr();
    values[0] *= RAD_TO_DEG;
    values[1] *= RAD_TO_DEG;
    assertVectorEqual(values, coords_actual, size, __LINE__);
}

void test_Coordinates_transformTo_update_wrapped_state(){
	// This is the standard PROJ.4 definition of a "WGS84" coordinate system.
	CoordinateSystem coordsys("+proj=longlat +ellps=WGS84 +vunits=m");
	
	int size = 6;
	// These are the coordinate values: {x0, y0, z0, x1, y1, z1}. The ordering of the x/y/z does not 
	// matter within a set but should be consistent across sets.
	double values[] = {290, 40, 0, 360, 50, 1};

	// If there is no z-coordinate, use -1 instead of the index location (2 in this case).
	TestData td(values, size, 0, 1, 2, coordsys);
	
	// This is the destination coordinate system. The default definition for a spherical coordinate
	// system. Note the conversion of vertical units to kilometers.
	CoordinateSystem coordsys_dest("+proj=longlat +a=6370997 +b=6370997 +vunits=km");
	
	assertTrue(coordsys_dest.getWrappedState() == UNINITIALIZED);
	// This transforms the coordinate values to match the destination coordinate system.
	td.p_coords->transformTo(coordsys_dest);
	assertTrue(coordsys_dest.getWrappedState() == coordsys.getWrappedState());
	
	// Convert back to degrees remembering to skip the z-coordinate.
    double *p_values = td.p_coords->getValuesPtr();
    for (int i = 0; i < size; i++){
		if (i != 2 && i != 5) {
			p_values[i] *= RAD_TO_DEG;
		}
	}
	
	double values_actual[] = {290, 39.8106, -1.65214, 360, 49.8104, -5.36448};
	assertVectorEqual(p_values, values_actual, size, __LINE__);
}

void test_proj4StringToMap_proj4MapToString(){
    int size = 2;
    vector<string> sources(size, "");
    sources[0] = "+proj=longlat +ellps=wgs84 +no_defs ";
    sources[1] = "+proj=longlat +ellps=wgs84 +no_defs";
    proj_map actual;
    proj_map pm;
    actual["proj"] = "longlat";
    actual["ellps"] = "wgs84";
    actual["no_defs"] = "";
    bool cmp;
    string proj_string;

    for (int i = 0; i < size; i++){
        pm = proj4StringToMap(sources[i]);
        cmp = pm.size() == actual.size() && equal(pm.begin(), pm.end(), actual.begin());
        assertTrue(cmp, __LINE__);
        proj_string = proj4MapToString(pm);
        cmp = proj_string == "+ellps=wgs84 +no_defs +proj=longlat ";
        assertTrue(cmp, __LINE__);
    }
}

void runTest(void (*f)(), string name){
    TEST_TOTAL += 1;
    try{
        f();
        cout << "%TEST_SUCCESS% " << name << " (" << __FILE__ << ")" << endl;
    }
    catch (SkipTest& e){
        TEST_SKIPPED += 1;
        cout << "%TEST_SKIPPED% " << name << " (" << __FILE__ << ")" << endl;
        cout << e.what();
    }
    catch (exception& e){
        TEST_FAILED += 1;
        cout << "%TEST_FAILED% " << name << " (" << __FILE__ << ")" << endl;
        cout << e.what();
    }
}

int main(int argc, char** argv) {
    TEST_TOTAL = 0;
    TEST_SKIPPED = 0;
    TEST_FAILED = 0;

    cout << "%SUITE_STARTED% " << __FILE__ << endl;

    runTest(&test_Coordinates, "test_Coordinates");
    runTest(&test_Coordinates_updateWrappedState, "test_Coordinates_updateWrappedState");
    runTest(&test_Coordinates_transformTo, "test_Coordinates_transformTo");
    runTest(&test_Coordinates_transformTo_bigger_arrays, "test_Coordinates_transformTo_bigger_arrays");
    runTest(&test_Coordinates_transformTo_edge_effects, "test_Coordinates_transformTo_edge_effects");
    runTest(&test_Coordinates_transformTo_edge_effects_wrapped, "test_Coordinates_transformTo_edge_effects_wrapped");
    runTest(&test_Coordinates_transformTo_rotated_pole, "test_Coordinates_transformTo_rotated_pole");
    runTest(&test_Coordinates_transformTo_three_dimensions, "test_Coordinates_transformTo_three_dimensions");
    runTest(&test_Coordinates_transformTo_update_wrapped_state, "test_Coordinates_transformTo_update_wrapped_state");
    runTest(&test_Coordinates_transformTo_various_transformations,
            "test_Coordinates_transformTo_various_transformations");
    runTest(&test_CoordinateSystem, "test_CoordinateSystem");
    runTest(&test_getWrappedState, "test_getWrappedState");
    runTest(&test_proj4StringToMap_proj4MapToString, "test_proj4StringToMap_proj4MapToString");

    cout << "%SUITE_FINISHED%" << endl << endl;
    cout << "  Total = " << TEST_TOTAL << endl;
    cout << "Skipped = " << TEST_SKIPPED << endl;
    cout << " Failed = " << TEST_FAILED << endl;

    int ret = EXIT_SUCCESS;
    if (TEST_FAILED != 0){
        ret = EXIT_FAILURE;
    }
    return (ret);
}
