#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include "transform.h"

using namespace std;

//todo consider inlining functions
//todo remove as many functions from transform.h as possible - like proj4MapContainsAll

int applyLongLatInitMap(const coord_overloads& cmap, Coordinates& coords){
    //todo doc
    map<int, double>::const_iterator it;
    double* values = coords.getValuesPtr();
    for (it = cmap.begin(); it != cmap.end(); it++){
        values[it->first] = it->second;
    }
    return 0;
};

coord_overloads getLongLatInitMap(const Coordinates& coords){
	//todo doc
    coord_overloads cmap;
    double edge = 180*DEG_TO_RAD;
	const CoordinateIndex* cindex = coords.getCoordinateIndexPtr();
    double* values = coords.getValuesPtr();
	
    for (int i = cindex->getXLocation(); i < coords.getSize(); i+= coords.getNdim()){
        double coord_current = values[i];
        if (coord_current == 0){
            cmap[i] = 0;
        }
        else if (abs(coord_current) == edge){
            cmap[i] = coord_current;
        }
    }
    return cmap;
};

WrappedState getWrappedState(const Coordinates& coords){
    // TODO doc
    WrappedState result = UNDEFINED;
    double coord_current;
	const CoordinateIndex* cindex = coords.getCoordinateIndexPtr();
    double* values = coords.getValuesPtr();

    for (int i = cindex->getXLocation(); i < coords.getSize(); i+=coords.getNdim()){
        coord_current = values[i];
        if (coord_current > 180*DEG_TO_RAD){
            result = UNWRAPPED;
            break;
        }
        if (coord_current < 0){
            result = WRAPPED;
            break;
        }
    }

    return result;
};

bool proj4MapContainsAll(const proj_map& pm, const string& proj_name, const vector<string>& keys){
    bool result = true;

    if (pm.find("proj")->second != proj_name){
        result = false;
    }
    else{
        for (vector<string>::const_iterator it = keys.begin(); it != keys.end(); ++it){
            if (pm.find(*it) == pm.end()) {
                result = false;
                break;
            }
        }
    }

    return result;
};

string proj4MapToString(const proj_map& pm){
    string result;
    for (map<string, string>::const_iterator it = pm.begin(); it != pm.end(); it++){
        result += "+" + it->first;
        if (it->second != ""){
            result += "=" + it->second;
        }
        result += " ";
    }
    return result;
};

proj_map proj4StringToMap(const string& source){
    proj_map result;
    vector<size_t> locations;
    string delimiter = " ";
    vector<string> tokens;

    int start = 0;
    int end = 0;
    while ((end = source.find(delimiter, start)) != string::npos) {
        locations.push_back(end);
        start = end + 1;
    }

    int pos = 0;
    for (vector<size_t>::iterator it = locations.begin(); it != locations.end(); ++it) {
        tokens.push_back(source.substr(pos, *it - pos));
        pos = *it + 1;
    }
    tokens.push_back(source.substr(pos));

    string key;
    string value;
    string current;
    size_t parm_delimiter;
    size_t value_delimiter;
    for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it){
        current = *it;
        parm_delimiter = current.find("+");
        if (parm_delimiter != string::npos){
            value_delimiter = current.find("=");
            if (value_delimiter != string::npos){
                key = current.substr(1, value_delimiter - 1);
                value = current.substr(value_delimiter + 1);
            }
            else{
                key = current.substr(1, string::npos);
                value = "";
            }
            result[key] = value;
        }
    }

    return result;
};

bool stringContainsAll(const string& target, const vector<string>& strings){
    bool result;
    vector<bool> contains(strings.size(), true);

    for(vector<string>::size_type i = 0; i != strings.size(); i++) {
        size_t has_string = target.find(strings[i]);
        if (has_string == string::npos){
            contains[i] = false;
        }
    }

    result = true;
    for(vector<bool>::iterator it = contains.begin(); it != contains.end(); ++it) {
        if (!*it){
            result = false;
            break;
        }
    }

    return result;
};

// #####################################################################################################################

CoordinateSystem::CoordinateSystem(const string& proj_string){
    this->initialize(proj_string, UNINITIALIZED, true);
};

CoordinateSystem::CoordinateSystem(const string& proj_string, WrappedState wrapped_state){
    this->initialize(proj_string, wrapped_state, true);
};

CoordinateSystem::CoordinateSystem(const string& proj_string, WrappedState wrapped_state, bool should_update_validate){
    this->initialize(proj_string, wrapped_state, should_update_validate);
};

//todo add initialization with coordinate system map for CoordinateSystem

CoordinateSystem::~CoordinateSystem() {
    pj_free(this->proj);
};

string CoordinateSystem::asString() const {
	return proj4MapToString(this->value);
};

bool CoordinateSystem::containsToWGS84() const {
    bool result = true;
    if (this->value.find("towgs84") == this->value.end()) {
        result = false;
    }
    return result;
};

void CoordinateSystem::initialize(const string& proj_string, WrappedState wrapped_state, bool should_update_validate){
    this->value = proj4StringToMap(proj_string);
    this->wrapped_state = wrapped_state;

    // Modify the internal proj4 map to include some default parameters.
    if (should_update_validate){
        // http://lists.osgeo.org/pipermail/mapserver-users/2003-November/046863.html
        if (this->value.find("no_defs") == this->value.end()) {
            this->value["no_defs"] = "";
        }

        if (this->isSphericalLongLat() || this->isWGS84LongLat()){
			vector<string> keys(1, "lon_wrap");
			bool contains_lonwrap = proj4MapContainsAll(this->value, "longlat", keys);
			
            if (!this->containsToWGS84()){
                this->value["towgs84"] = ESMF_PROJ4_TOWGS84;
            }
			if (!contains_lonwrap){
				if (this->wrapped_state == UNWRAPPED){
					this->value["lon_wrap"] = "180";
				}
			}
			else {
				if (this->value["lon_wrap"] == "180"){
					this->wrapped_state = UNWRAPPED;
				}
			}
        }
    }

    const char* cproj_string = proj4MapToString(this->value).c_str();
    this->proj = pj_init_plus(cproj_string);
};

bool CoordinateSystem::isLongLat() const {
    return pj_is_latlong(this->proj);
};

bool CoordinateSystem::isSphericalLongLat() const {
    vector<string> keys(2, "");
    keys[0] = "a";
    keys[1] = "b";

    return proj4MapContainsAll(this->value, "longlat", keys);
};

bool CoordinateSystem::isWGS84LongLat() const {
    // todo add map comparison that uses a subset as opposed to full equality
    // todo what to do about the capital WGS84? if it is not capital, a segfault will occur.
    bool result;
    string proj_string = proj4MapToString(this->value);

    vector<string> strings1(2, "");
    strings1[0] = "+proj=longlat";
    strings1[1] = "+ellps=WGS84";

    vector<string> strings2(2, "");
    strings2[0] = "+proj=longlat";
    strings2[1] = "+datum=WGS84";

    bool has1 = stringContainsAll(proj_string, strings1);
    bool has2 = stringContainsAll(proj_string, strings2);

    if (has1 || has2){
        result = true;
    }
    else{
        result = false;
    }

    return result;
};

bool CoordinateSystem::isWrappableLongLat() const {
    bool result = false;
    if (this->isLongLat() && (this->isSphericalLongLat() || this->isWGS84LongLat())){
        result = true;
    }
    return result;
};

void CoordinateSystem::setWrappedState(WrappedState wrapped_state){
    pj_free(this->proj);
    this->initialize(proj4MapToString(this->value), wrapped_state, true);
};

bool CoordinateSystem::operator==(const CoordinateSystem& rhs) const {
    const proj_map& lhs_value = this->value;
    const proj_map& rhs_value = rhs.value;
    return lhs_value.size() == rhs_value.size() && equal(lhs_value.begin(), lhs_value.end(), rhs_value.begin());
};

bool CoordinateSystem::operator!=(const CoordinateSystem& rhs) const {
    return !(this == &rhs);
};

// #####################################################################################################################

int CoordinateIndex::getDimensionCount() const {
    int result;

    if (this->z == -1){
        result = 2;
    }
    else{
        result = 3;
    }
    return result;
};

// #####################################################################################################################

Coordinates::Coordinates(double values[], int size, CoordinateIndex& coord_index, CoordinateSystem& coordsys) {
    this->p_values = values;
    this->size = size;
    this->p_coordsys = &coordsys;
    this->p_coord_index = &coord_index;
    this->ndim = coord_index.getDimensionCount();
    
	if (coordsys.getWrappedState() == UNINITIALIZED && coordsys.isLongLat()){
		this->updateWrappedState();
	}
};

void Coordinates::updateWrappedState(){
    WrappedState wrapped_state = getWrappedState(*this);
    this->p_coordsys->setWrappedState(wrapped_state);
};

int Coordinates::transformTo(CoordinateSystem& coordsys_dest){
	//todo add option to not get/apply init maps
    int p;
    double *p_z;
    p_z = NULL;
	coord_overloads cmap;
	
	bool is_longlat = this->p_coordsys->isLongLat() && coordsys_dest.isLongLat();
	
	if (is_longlat){
		if (coordsys_dest.getWrappedState() == UNINITIALIZED){
			coordsys_dest.setWrappedState(this->p_coordsys->getWrappedState());
		}
		cmap = getLongLatInitMap(*this);
	}

    // Point to the third dimension (z) coordinate.
    if (this->ndim != 2){
        p_z = &this->p_values[this->p_coord_index->getZLocation()];
    }

    p = pj_transform(*this->p_coordsys->getProjPJPtr(), *coordsys_dest.getProjPJPtr(), this->size/this->ndim, this->ndim,
                     &this->p_values[this->p_coord_index->getXLocation()], &this->p_values[this->p_coord_index->getYLocation()], 
					 p_z);

    this->p_coordsys = &coordsys_dest;

	if (is_longlat){
		applyLongLatInitMap(cmap, *this);
	}

    return p;
};
