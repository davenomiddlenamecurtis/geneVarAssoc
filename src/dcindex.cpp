#if 0
Copyright 2018 David Curtis

This file is part of the geneVarAssoc package.

geneVarAssoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

geneVarAssoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with geneVarAssoc.If not, see <http://www.gnu.org/licenses/>.
#endif

/* this index will return 0L if it does not find a matching record */
 
#include "dcerror.hpp"
#include "dcindex.hpp"

int dcIndex::add(char *key,long rec)
// not checking here to see if already exists
{
	m.insert(std::pair<std::string, long>(key, rec));
	return 1;
}

long dcIndex::get_last()
{
	return (m.empty()) ? 0L : (it= --m.end())->second;
}


int dcIndex::remove()
{
	if (m.empty())
		return 0;
	else
	{
		it=m.erase(it);
		return 1;
	}
}

int dcIndex::remove(char *key)
{
	return m.erase(key);
}

long dcIndex::get_next()
{
	++it;
	return (it == m.end()) ? 0L : it->second;
}

long dcIndex::get_prev()
{
	if (m.empty() || it == m.begin())
		return 0L;
	else
		return --it->second;
}

long dcIndex::get_first()
{
	if (m.empty())
		return 0L;
	else 
	{
		it = m.begin();
		return it->second;
	}
}

int dcIndex::find_matching_node(char *key,long rec)
{
// find first exact match, then try going backwards and forwards from it
// till hit right record number
	if (m.empty() || (it=m.find(key))==m.end())
		return 0L;
	if (it->second != rec)
	{
		dcerror(1, "error in int dcIndex::find_matching_node(char *key,long rec): key=%s, rec=%ld but there is an entry for that key with rec=%ld\n",
			key, rec, it->second);
		return 0L;
	}
}


long dcIndex::exact_find(char *key)
{
	if (m.empty() || (it = m.find(key)) == m.end())
		return 0L;
	else
		return it->second;
}

long dcIndex::near_find(char *key)
{
	if (m.empty() )
		return 0L;
	else
	{
		if ((it = m.lower_bound(key)) == m.end())
			it = m.begin();
		return it->second;
	}
}

int dcIndex::make_new(char *name)
{
	fn = name;
	m.clear();
	return 1;
}

int dcIndex::open_old(char *name)
{
	FILE *fp;
	char buff[MAXKEYLENGTH+1],*ptr;
	long rec;
	fp = fopen(name, "r");
	if (fp==0)
	{
		dcerror(1, "Could not open index file %s", name);
		return 0;
	}
	m.clear();
	fn = name;
	while (1)
	{
		if (!fgets(buff, MAXKEYLENGTH, fp))
			break;
		ptr = strchr(buff, '\n');
		if (ptr == 0)
		{
			fclose(fp);
			dcerror(1, "This line is too long in index file %s:\n%s\n", name, buff);
			return 0;
		}
		else 
			*ptr = '\0';
		if (!fgets(buff, MAXKEYLENGTH, fp) || (rec=atol(buff),rec==0))
		{
			fclose(fp);
			dcerror(1, "Could not read valid record number from index file %s for this key:\n%s\n", name, buff);
			return 0;
		}
		m.insert(std::pair<std::string, long>(buff, rec));
	}
	fclose(fp);
	it = m.begin();
	return 1;
}

dcIndex::~dcIndex()
{
	close();
}

void dcIndex::close()
{
	FILE *fp;
	fp = fopen(fn.c_str(), "w");
	if (fp == 0)
		dcerror(1, "Could not open index file %s", fn.c_str());
	else
	{
		for (it = m.begin(); it != m.end(); ++it)
			fprintf(fp, "%s\n%ld\n", it->first.c_str(), it->second);
	}
	fclose(fp);
	m.clear();
	fn = "";
}

dcIndex::dcIndex()
{ 
	fn = "";
}

