#ifndef __TRANSFORM_H_INCLUDED__
#define __TRANSFORM_H_INCLUDED__

// #####################################################################################################################

#include <proj_api.h>

#include <map>
#include <string>
#include <vector>

// #####################################################################################################################

// State indicators for the wrapping of a longitude/latitude coordinate system. Wrapped means the longitude
// values occur in -180 to 180 domain. Unwrapped indicates 0 to 360. Undefined indicates 0 to 180. Unintialized
// means the wrapping has not been set.
enum WrappedState {WRAPPED=0, UNWRAPPED=1, UNDEFINED=2, UNINITIALIZED=3};
// Holds coordinate values to replace to manage PROJ.4 edge effects. For example, in an unwrapped coordinate
// system, 0 is remapped to 360 by PROJ.4.
typedef std::map<int, double> coord_overloads;
// Holds key-value representations of PROJ.4 parameter-values.
typedef std::map<std::string, std::string> proj_map;
// The default "towgs84" parameter values for PROJ.4 definitions.
const std::string ESMF_PROJ4_TOWGS84 = "0,0,0";

// #####################################################################################################################

class Uncopyable {
    protected:
        Uncopyable(){};
        ~Uncopyable(){};
    private:
        Uncopyable(const Uncopyable&){};
        Uncopyable& operator=(const Uncopyable&){};
};

// TODO make data members private
class CoordinateSystem: private Uncopyable {
    public:
        CoordinateSystem(const std::string&);
        CoordinateSystem(const std::string&, WrappedState);
        CoordinateSystem(const std::string&, WrappedState, bool);

		std::string asString() const;
        bool containsToWGS84() const;
		const projPJ* getProjPJPtr() const {
			return &this->proj;
		}
		const proj_map* getValuePtr() const {
			return &this->value;
		}
		WrappedState getWrappedState() const {
			return this->wrapped_state;
		}
        bool isLongLat() const;
        bool isSphericalLongLat() const;
        bool isWGS84LongLat() const;
        bool isWrappableLongLat() const;
        void setWrappedState(WrappedState);

		bool operator==(const CoordinateSystem&) const;
        bool operator!=(const CoordinateSystem&) const;
		
        ~CoordinateSystem();
    private:
		proj_map value;
		projPJ proj;
		WrappedState wrapped_state;
		
        void initialize(const std::string&, WrappedState, bool);
};

class CoordinateIndex: private Uncopyable {
	public:
		CoordinateIndex(int x, int y, int z){
			this->x = x;
			this->y = y;
			this->z = z;
		};
		
		int getDimensionCount() const;
		int getXLocation() const {return this->x;}
		int getYLocation() const {return this->y;}
		int getZLocation() const {return this->z;}
	private:
		int x;
		int y;
		int z;
};

class Coordinates: private Uncopyable {
public:
        Coordinates(double values[], int size, CoordinateIndex& coord_index, CoordinateSystem& coordsys);

		const CoordinateIndex* getCoordinateIndexPtr() const {return this->p_coord_index;}
		const CoordinateSystem* getCoordinateSystemPtr() const {return this->p_coordsys;}
		int getNdim() const {return this->ndim;}
		int getSize() const {return this->size;}
		double * getValuesPtr() const {return this->p_values;}
        void updateWrappedState();
        int transformTo(CoordinateSystem&);
	private:
		CoordinateIndex* p_coord_index;
		double* p_values;
		int size;
		int ndim;
		CoordinateSystem* p_coordsys;
};

// #####################################################################################################################

int applyLongLatInitMap(const coord_overloads& cmap, Coordinates& coords);
coord_overloads getLongLatInitMap(const Coordinates& coords);
WrappedState getWrappedState(const Coordinates& coords);
bool proj4MapContainsAll(const proj_map& pm, const std::string& proj_name, const std::vector<std::string>& keys);
std::string proj4MapToString(const proj_map& pm);
proj_map proj4StringToMap(const std::string& source);
bool stringContainsAll(const std::string& target, const std::vector<std::string>& strings);

// #####################################################################################################################

#endif
