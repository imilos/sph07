#ifndef MMATERIALDATA_H_
#define MMATERIALDATA_H_

#include "global.h"
#include "MMaterial.h"
#include <vector>

class MMaterialData
{
public:
	MMaterialData();
	virtual ~MMaterialData();
	
public:
	int AddMaterial(const MMaterial &mat);
	int RemoveMaterial(int ID);
	MMaterial& GetMaterial(int ID);
	int DeleteContents();
	int GetSize();
	MMaterial& operator[](int ID);
	
private:
	std::vector<MMaterial> m_MatList;
};

#endif /*MMATERIALDATA_H_*/
