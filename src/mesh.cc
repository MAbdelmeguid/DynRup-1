#include "mesh.h"



void Mesh::Abaqus_IO(string &fname)
{
    std::string s;
        std::ifstream _in(fname);
        int mm = 0;
        while (true)
        {
            std::getline(_in,s);
            // Convert s to uppercase
            std::string upper(s);
            std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
            // 0.) Look for the "*Part" Section
    //        if (upper.find("*PART")== static_cast<std::string::size_type>(0))
    //        {
    //            std::cout<<"Find *PART"<<std::endl;
    //        }
            // 1.) Loo for the "*Nodes" section
            if (upper.find("*NODE")== static_cast<std::string::size_type>(0))
            {
                //std::string nset_name = s ;
                //std::cout<<s<<std::endl;
                // Temperatry variables for parsing lines of text
                char c ;
                std::string line;
                int i=0;
                while (_in.peek()!='*'&&_in.peek()!=EOF)
                {
                    // Read an entire line which corresponds to a single points's id and (x,y) value
                    std::getline(_in,line);
                    // Revomie all whitesspaces characters from the line
                    line.erase(std::remove_if(line.begin(),line.end(),::isspace),line.end());
                    // Make a stream out of the modified line so we can stream values from it in the usaly way
                    std::stringstream ss(line);
                    int abaqus_node_id = 0;
                    double x = 0 , y = 0;
                    ss >> abaqus_node_id >> c >> x >>c >> y ;
                    Node.conservativeResize(i+1, 2);
                    Node.row(i)<< x , y ;
                    i++;
                }
                //  std::cout<< Nodes<< "\n"<< Nodes.rows()<<std::endl;
            }
            else if (upper.find("*ELEMENT,")==static_cast<std::string::size_type>(0))
            {
                //std::string elset_name = s;
                //std::cout<<s <<std::endl;
                char c ;
                std::string line;
                int i = 0;
                while (_in.peek()!='*'&&_in.peek()!=EOF)
                {
                    // Read an entire line which corresponds to a single element's id and connectivity value (Q4_only)
                    std::getline(_in,line);
                    // Revomie all whitesspaces characters from the line
                    line.erase(std::remove_if(line.begin(),line.end(),::isspace),line.end());
                    // Make a stream out of the modified line so we can stream values from it in the usaly way
                    std::stringstream ss(line);
                    int abaqus_el_id = 0;
                    int node1 =0 , node2 =0 , node3 = 0 , node4 =0;
                    ss >> abaqus_el_id >> c >> node1 >> c >>node2 >> c >> node3 >> c >> node4 ;
                    // add -1 here becasue we are using a zero node numbering zero is the first node
                    Element.conservativeResize(i+1, 4);
                    Element.row(i)<<node1-1,node2-1,node3-1,node4-1;
                    i++;
                }
            }
            if (_in.eof())
                break;
        }
    }
