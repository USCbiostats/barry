#include<regex>
#include<iostream>
#include<string>
#include<iterator>
#include<vector>

int main() {

    size_t m_order = 1;
    size_t y_ncol  = 4;

    std::string formula1("{y1, y3}");

    /**
     * Matrix representation
     * |n/a 1 n/a 1|
     * |n/a 1 n/a 0|
     * 
     * Vector
     * {2, 3, 6, 7}
     * {true, true, true, false}
     */
    std::string formula2("{y1_0, y3} > {y1_1, 0y3_1}");
    std::string formula3("{y1_0}");

    std::regex model_cross("\\{\\s*0?y[0-9]+(\\s*,\\s*0?y[0-9]+)*\\s*\\}");
    std::regex model_long("\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\}\\s*(>)\\s*\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\s*\\}");

    auto empty = std::sregex_iterator();

    std::vector< size_t > res_locations;
    std::vector< size_t > res_sign;

    // This column-major vector indicates true if the variable has already been
    // selected
    std::vector< bool > selected((m_order + 1) * y_ncol, false);

    for (auto s : {formula1, formula2, formula3})
    {

        std::smatch match;
        std::regex_search(s, match, model_long);
        if (!match.empty())
        {

            // Will indicate where the arrow is located at
            size_t arrow_position = match.position(4u);

            // This pattern will match 
            std::regex pattern;
            if (m_order > 1)
                pattern = ("(0?)y([0-9]+)_([0-9]+)");
            else
                pattern = ("(0?)y([0-9]+)(_[0-9]+)?");

            auto iter = std::sregex_iterator(s.begin(), s.end(), pattern);

            for (auto i = iter; i != empty; ++i)
            {

                // Baseline position
                size_t current_location = i->position(0u);

                // First value true/false
                bool is_positive;
                if (i->operator[](1u).str() == "")
                    is_positive = true;
                else if (i->operator[](1u).str() == "0")
                    is_positive = false;
                else
                    throw std::logic_error("The number preceding y should be either none or zero.");

                // Variable position
                size_t y_col = std::stoul(i->operator[](2u).str());
                if (y_col >= y_ncol)
                    throw std::logic_error("The proposed column is out of range.");

                // Time location
                size_t y_row;
                if (m_order > 1)
                {

                    y_row = std::stoul(i->operator[](3u).str());

                    if (y_row >= m_order)
                        throw std::logic_error("The proposed row is out of range.");

                } else
                    y_row = (i->position() < arrow_position ? 0u: 1u);

                if (selected[y_col * (m_order + 1) + y_row])
                    throw std::logic_error(
                        "The term " + i->str() + " shows more than once in the formula.");

                // Only the end of the chain can be located at position after the
                // arrow
                if ((current_location > arrow_position) && (y_row != m_order))
                    throw std::logic_error(
                        "Only the row " + std::to_string(m_order) +
                        " can be specified at the RHS of the motif."
                        );


                selected[y_col * (m_order + 1) + y_row] = true;

                res_locations.push_back(y_col * (m_order + 1) + y_row);
                res_sign.push_back(is_positive);
                

            }

            auto b = std::sregex_iterator(s.begin(), s.end(), pattern);
            for (auto i = b; i != std::sregex_iterator(); ++i)
            {
                std::smatch m_match = *i;
                std::cout << m_match.str() << std::endl;
                for (auto j = m_match.begin(); j != m_match.end(); ++j)
                    std::cout<< " - " << j->str() << " at " << m_match.position() << std::endl;
            }


        } else if (std::regex_match(s, model_cross)){

            std::cout << "is cross-sectional" << std::endl;

        } else 
            std::cout << "Syntax error in DEFM" << std::endl;
    }


}