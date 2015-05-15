#include <transform.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv){	
	// This is the standard PROJ.4 definition of a "WGS84" coordinate system. When the code to read coordinate systems
    // directory from netCDF attributes is complete, a CoordinateSystem object will be returned from the parser negating
    // the need to work with strings.
	CoordinateSystem coordsys("+proj=longlat +ellps=WGS84 +vunits=m");
	// Some parameters are added by the software. Ordering of parameters do not matter in the PROJ.4 string.
	cout << coordsys.asString() << endl;
	
	int size = 6;
	// These are the coordinate values ordered by: {x0, y0, z0, x1, y1, z1}. The ordering of the x/y/z does not matter
	// within a set but should be consistent across sets. A set consists of a coordinate pair or triple.
	double values[] = {290, 40, 0, 360, 50, 1};
	// Convert the coordinate values to radians skipping the z-axis.
	for (int i = 0; i < size; i++){
		if (i != 2 && i != 5) {
            // This is the conversion value built-in to PROJ.4.
			values[i] *= DEG_TO_RAD;
		}
	}
	
	// If there is no z-coordinate, use -1 instead of the index location (2 in this case). This object may be extraneous
    // but it reduces the argument count to Coordinates.
	CoordinateIndex cindex(0, 1, 2);
	Coordinates coords(values, size, cindex, coordsys);
	
	// This is the destination coordinate system. The default definition for a spherical coordinate system. Note the
	// conversion of vertical units to kilometers.
	CoordinateSystem coordsys_dest("+proj=longlat +a=6370997 +b=6370997 +vunits=km");
    cout << coordsys_dest.asString() << endl;
	
	// This transforms the coordinate values to match the destination coordinate system.
	coords.transformTo(coordsys_dest);
	
	cout << endl << "----- Transformed Coordinates:" << endl;
	// Convert back to degrees remembering to skip the z-coordinate.
	double* p_values = coords.getValuesPtr();
	for (int i = 0; i < size; i++){
		if (i != 2 && i != 5) {
			p_values[i] *= RAD_TO_DEG;
		}
		cout << p_values[i] << endl;
	}

    return EXIT_SUCCESS;
};