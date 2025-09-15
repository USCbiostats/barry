#ifndef BARRY_DEFM_MOTIF_FORMULA_HPP
#define BARRY_DEFM_MOTIF_FORMULA_HPP
/**
 * @brief Parses a motif formula
 * 
 * @details This function will take the formula and generate the corresponding
 * input for defm::counter_transition(). Formulas can be specified in the
 * following ways:
 * 
 * - Intercept effect: {...} No transition, only including the current state.
 * - Transition effect: {...} > {...} Includes current and previous states.
 * 
 * The general notation is `[0]y[column id]_[row id]`. A preceeding zero
 * means that the value of the cell is considered to be zero. The column
 * id goes between 0 and the number of columns in the array - 1 (so it
 * is indexed from 0,) and the row id goes from 0 to m_order.
 * 
 * ## Intercept effects
 * 
 * Intercept effects only involve a single set of curly brackets. Using the
 * 'greater-than' symbol (i.e., '>') is only for transition effects. When
 * specifying intercept effects, users can skip the `row_id`, e.g.,
 * `y0_0` is equivalent to `y0`. If the passed `row id` is different from
 * the Markov order, i.e., `row_id != m_order`, then the function returns
 * with an error. 
 * 
 * Examples:
 * 
 * - `"{y0, 0y1}"` is equivalent to set a motif with the first element equal
 * to one and the second to zero. 
 * 
 * ## Transition effects
 * 
 * Transition effects can be specified using two sets of curly brackets and
 * an greater-than symbol, i.e., `{...} > {...}`. The first set of brackets,
 * which we call LHS, can only hold `row id` that are less than `m_order`.
 * 
 * 
 * @param formula A string specifying the motif formula (see details).
 * @param locations A vector of locations for the motif variables.
 * @param signs A vector of signs for the motif variables.
 * @param m_order The Markov order.
 * @param y_ncol The number of columns in the response variable.
 * @param covar_name A string to hold the name of the covariate (if any).
 * @param vname A string to hold the variable name (if any).
 */
inline void defm_motif_parser(
    std::string formula,
    std::vector< size_t > & locations,
    std::vector< bool > & signs,
    size_t m_order,
    size_t y_ncol,
    std::string & covar_name,
    std::string & vname
)
{
    // Resetting the results
    locations.clear();
    signs.clear();

    std::regex pattern_intercept(
        std::string("\\{\\s*[01]?y[0-9]+(_[0-9]+)?(\\s*,\\s*[01]?y[0-9]+(_[0-9]+)?)*\\s*\\}") +
        std::string("(\\s*x\\s*[^\\s]+([(].+[)])?\\s*)?")
        );
    std::regex pattern_transition(
        std::string("\\{\\s*[01]?y[0-9]+(_[0-9]+)?(\\s*,\\s*[01]?y[0-9]+(_[0-9]+)?)*\\}\\s*(>)\\s*") +
        std::string("\\{\\s*[01]?y[0-9]+(_[0-9]+)?(\\s*,\\s*[01]?y[0-9]+(_[0-9]+)?)*\\s*\\}") +
        std::string("(\\s*x\\s*[^\\s]+([(].+[)])?\\s*)?")
        );

    auto empty = std::sregex_iterator();

    // This column-major vector indicates true if the variable has already been
    // selected
    std::vector< bool > selected((m_order + 1) * y_ncol, false);

    std::smatch match;
    std::regex_match(formula, match, pattern_transition);
    if (!match.empty())
    {

        if (m_order == 0)
            throw std::logic_error("Transition effects are only valid when the data is a markov process.");

        // Matching the pattern '| [no spaces]$'
        std::regex pattern_conditional(".+[}]\\s+x\\s+([^(]+)([(][^)]+[)])?\\s*$");
        std::smatch condmatch;
        std::regex_match(formula, condmatch, pattern_conditional);
        // Extracting the [no_spaces] part of the conditional
        if (!condmatch.empty())
        {
            covar_name = condmatch[1].str();
            vname = condmatch[2].str();

            // Removing starting and ending parenthesis
            if (vname != "")
                vname = vname.substr(1, vname.size() - 2);

        }

        // Will indicate where the arrow is located at
        size_t arrow_position = match.position(4u);

        // This pattern will match 
        std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

        auto iter = std::sregex_iterator(formula.begin(), formula.end(), pattern);

        for (auto i = iter; i != empty; ++i)
        {

            // Baseline position
            size_t current_location = i->position(0u);

            // First value true/false
            bool is_positive = true;
            if (i->operator[](1u).str() == "0")
                is_positive = false;

            // Variable position
            size_t y_col = std::stoul(i->operator[](2u).str());
            if (y_col >= y_ncol)
                throw std::logic_error("The proposed column is out of range.");

            // Time location
            size_t y_row;
            std::string tmp_str = i->operator[](4u).str();
            if (m_order > 1)
            {
                // If missing, we replace with the location 
                if (tmp_str == "")
                {

                    if (current_location > arrow_position)
                        y_row = m_order;
                    else
                        throw std::logic_error("LHS of transition must specify time when m_order > 1");

                } else
                    y_row = std::stoul(tmp_str);

                if (y_row > m_order)
                    throw std::logic_error("The proposed row is out of range.");


            } else {

                // If missing, we replace with the location 
                if (tmp_str != "")
                    y_row = std::stoul(tmp_str);
                else
                    y_row = (current_location < arrow_position ? 0u: 1u);

            }

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

            locations.push_back(y_col * (m_order + 1) + y_row);
            signs.push_back(is_positive);
            

        }

        return;

    } 
    
    std::regex_match(formula, match, pattern_intercept);
    if (!match.empty())
    {

        // Matching the pattern '| [no spaces]$'
        std::regex pattern_conditional(".+[}]\\s+x\\s+([^(]+)([(][^)]+[)])?\\s*$");
        std::smatch condmatch;
        std::regex_match(formula, condmatch, pattern_conditional);
        // Extracting the [no_spaces] part of the conditional
        if (!condmatch.empty())
        {
            covar_name = condmatch[1].str();
            vname = condmatch[2].str();

            // Removing starting and ending parenthesis
            if (vname != "")
                vname = vname.substr(1, vname.size() - 2);
        }

        // This pattern will match 
        std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

        auto iter = std::sregex_iterator(formula.begin(), formula.end(), pattern);

        for (auto i = iter; i != empty; ++i)
        {
            
            // First value true/false
            bool is_positive = true;
            if (i->operator[](1u).str() == "0")
                is_positive = false;

            // Variable position
            size_t y_col = std::stoul(i->operator[](2u).str());
            if (y_col >= y_ncol)
                throw std::logic_error("The proposed column is out of range.");

            // Time location
            size_t y_row;
            if (i->operator[](4u).str() == "") // Assume is the last
                y_row = m_order;
            else {

                y_row = std::stoul(i->operator[](4u).str());

                if (y_row != m_order)
                    throw std::logic_error(
                        std::string("Intercept motifs cannot feature past events. ") +
                        std::string("Only transition motifs can: {...} > {...}.")
                        );

            }

            if (selected[y_col * (m_order + 1) + y_row])
                throw std::logic_error(
                    "The term " + i->str() + " shows more than once in the formula.");

            selected[y_col * (m_order + 1) + y_row] = true;

            locations.push_back(y_col * (m_order + 1) + y_row);
            signs.push_back(is_positive);
            

        }

        return;

    } 
    
    throw std::logic_error(
        "The motif specified in the formula: " + formula +
        " has the wrong syntax."
        );
    
}
#endif

