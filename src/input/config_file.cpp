// ==============================================================================
// This file is part of THOR.
//
//     THOR is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     THOR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     THOR directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Description: basic config file parser
//
//
// Method: 
//
// Known limitations: None.
//      
//
// Known issues: None.
//   
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// If you use this code please cite the following reference: 
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016  
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////
#include "config_file.h"

#include <regex>
#include <iomanip>
#include <fstream>
using namespace std;




config_file::config_file()
{
    append_config_var("config_version", version, -1);

}
        
bool config_file::parse_file(const string & filename)
{
    ifstream file(filename, std::ifstream::in);

    if (file.is_open())
    {
        return parse_config( file );
    }
    else
    {
        file.close();
        return false;
    }
}

    
bool config_file::parse_config(basic_istream<char> & config)
{
    
    // initialise key value pairs
    for (auto&&  it=config_vars.begin(); it!=config_vars.end(); ++it)
    {
        it->second->set_default();
        // cout << it->first  << "\t" << it->second->to_str() << endl;
    }
    
    


    // read input line by line
    regex comment_regex("^[[:blank:]]?#.*$");
    regex empty_line_regex("^[[:blank:]]*$");
    regex key_value_pair_regex("^[[:blank:]]*([A-Za-z0-9_]+)[[:blank:]]*="
                               "[[:blank:]]*(([[:alnum:]]|[[:punct:]])+)[[:blank:]]*(#.*)?$");
    
    int cnt = 1;
    for (std::string line; std::getline(config, line); ) {

        // cout << left << "line("<<setw(4)<<right<<cnt<<"): [" << line << "]" << setw(45 - line.length()) << " ";
        
        std::smatch match;
        
        if (regex_match(line, match, empty_line_regex))
        {
            // cout << "empty line" << endl;
        }
        else if (regex_match(line, match, comment_regex))
        {
            // cout << "comment"<<endl;
        }
        else if (regex_match(line, match, key_value_pair_regex))
        {
            // cout << "key value pair matched: (" << match[1] << ":" << match[2] <<")"<< endl;

            string key(match[1]);
            string value(match[2]);

            auto && it = config_vars.find(key);
            if (it != config_vars.end())
            {
                bool parsed = it->second->parse(value);
                if (!parsed)
                {
                     cout << "parsing of value ["<<value<<"] failed for key ["<<key<<"] "<<endl;                    
                }                
            }
            else
            {
                cout << "config file key ["<< key << "] does not exist" << endl;
                // return false;                
            }
        }
        else
        {
            cout << "error in config file at line " << cnt << endl;
            
            //return false;
        }
        cnt++;     
    }
    
    return true;
    
}
    
bool config_file::append_config_var(const string & name,
                                    std::unique_ptr<config_entry_interface> entry)
{
    // check if entry exists
    auto && it = config_vars.find(name);
    if (it != config_vars.end())
    {
        cout << "config entry " << name << " already exists"<< endl;
        return false;
    }
    else
    {
        // append entry
        config_vars[name] = std::move(entry);
        return true;
    }
    
}


