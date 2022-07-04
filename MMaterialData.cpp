#include "MMaterialData.h"

MMaterialData::MMaterialData()
{
}

MMaterialData::~MMaterialData()
{
}

int MMaterialData::AddMaterial(const MMaterial &mat)
{
	m_MatList.push_back(mat);
	return m_MatList.size()-1;
}

int MMaterialData::RemoveMaterial(int ID)
{
	// To be implemented
	return OK;
}

int MMaterialData::GetSize()
{
	return m_MatList.size();
}

MMaterial& MMaterialData::GetMaterial(int ID)
{
	return m_MatList[ID-1];
}

// IMPORTANT: The operator [] uses one based notation, not zero based
MMaterial& MMaterialData::operator[](int ID)
{
	return m_MatList[ID-1];
}

int MMaterialData::DeleteContents()
{
	m_MatList.clear();
	return OK;
}
