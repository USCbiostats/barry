#ifndef BARRY_DEFM_MOTIF_FORMULA_HPP
#define BARRY_DEFM_MOTIF_FORMULA_HPP
/**
 * @brief Parses a motif formula
 * 
 * @details This function will take the formula and generate the corresponding
 * input for defm::counter_generic(). Formulas can be specified in the
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
 * Transition effects can be specified using curly brackets separated by
 * greater-than symbols ('>').
 * 
 * **Two-group mode (backwards compatible):** `{...} > {...}`
 * - First group (LHS): Variables at times 0 to m_order-1. When m_order > 1,
 *   row indices must be explicitly specified.
 * - Second group (RHS): Variables at time m_order. Row indices can be omitted
 *   and will default to m_order.
 * 
 * **Multi-group mode (explicit):** `{...} > {...} > ... > {...}` (m_order+1 groups)
 * - Each group corresponds to a specific time point (0, 1, ..., m_order).
 * - Row indices can be omitted and will be inferred from group position.
 * - If specified, row indices must match the group position.
 * 
 * Examples:
 * - Order 1: `{y0_0} > {y0_1}` or `{y0_0} > {y0}` (both valid)
 * - Order 2: `{y0_0} > {y0}` (2-group, implicit final time)
 * - Order 2: `{y0_0} > {y0_1} > {y0_2}` (3-group, all explicit)
 *
 * ## Covariate interactions
 *
 * Both intercept and transition effects can be interacted with a covariate
 * by appending `x [covariate name]`, e.g., `{y0} > {y1} x Female`. The `x`
 * must be followed by whitespace. An optional label in parenthesis can
 * follow the covariate name. Motif variables are only read from the
 * bracketed groups, so covariate names such as `Day1` are allowed.
 *
 *
 * @param formula A string specifying the motif formula (see details).
 * @param locations A vector of locations for the motif variables.
 * @param signs A vector of signs for the motif variables.
 * @param covar_name If a covariate name is specified in the formula,
 * this variable will hold its name. If no covariate is specified, it will
 * be set to an empty string.
 * @param m_order The Markov order.
 * @param y_ncol The number of columns in the response variable.
 */
inline void defm_motif_parser(
    std::string formula,
    std::vector< size_t > & locations,
    std::vector< bool > & signs,
    size_t m_order,
    size_t y_ncol,
    std::string & covar_name
)
{
    // Resetting the results
    locations.clear();
    signs.clear();
    covar_name = "";

    // A single bracketed group, e.g., {y0, 0y1_1}. Only non-capturing groups
    // are used so the capture indices of the full patterns are predictable.
    const std::string pattern_group(
        "\\{\\s*[01]?y[0-9]+(?:_[0-9]+)?(?:\\s*,\\s*[01]?y[0-9]+(?:_[0-9]+)?)*\\s*\\}"
    );

    // Optional covariate interaction: 'x covar_name' followed by an optional
    // label in parenthesis. The 'x' must be followed by whitespace so that
    // the covariate name is unambiguous.
    // Captures: [2] = covariate name, [3] = label (with parenthesis).
    const std::string pattern_covar(
        "(?:\\s*x\\s+([^\\s(]+)\\s*([(].+[)])?)?\\s*"
    );

    // Captures: [1] = the motif part of the formula (the bracketed groups).
    std::regex pattern_intercept(
        "\\s*(" + pattern_group + ")" + pattern_covar
        );
    // Updated pattern to match one or more bracketed groups separated by '>'
    std::regex pattern_transition(
        "\\s*(" + pattern_group + "(?:\\s*>\\s*" + pattern_group + ")+)" +
        pattern_covar
        );

    auto empty = std::sregex_iterator();

    // This column-major vector indicates true if the variable has already been
    // selected
    std::vector< bool > selected((m_order + 1) * y_ncol, false);

    // Variables are only searched within the motif part of the formula, so
    // names in the covariate (e.g., 'Day1') are not taken as motif terms.
    std::string motif;
    std::string vname;
    auto extract_parts = [&](const std::smatch & m) -> void {

        motif      = m[1u].str();
        covar_name = m[2u].str();
        vname      = m[3u].str();

        // Removing starting and ending parenthesis
        if (vname != "")
            vname = vname.substr(1, vname.size() - 2);

    };

    std::smatch match;
    std::regex_match(formula, match, pattern_transition);
    if (!match.empty())
    {

        if (m_order == 0)
            throw std::logic_error("Transition effects are only valid when the data is a markov process.");

        extract_parts(match);

        // Find all bracketed groups to determine which time point each variable belongs to
        std::regex bracket_pattern("\\{[^}]+\\}");
        std::vector<std::pair<size_t, size_t>> bracket_ranges; // start, end positions

        auto brackets_begin = std::sregex_iterator(motif.begin(), motif.end(), bracket_pattern);
        for (auto i = brackets_begin; i != empty; ++i)
            bracket_ranges.push_back({i->position(), i->position() + i->length()});

        size_t num_groups = bracket_ranges.size();

        // For backwards compatibility, allow 2 groups (original behavior) or m_order+1 groups (new behavior)
        if (num_groups == 2)
        {
            // Two-group mode (backwards compatible):
            // - First group: variables at time 0 to m_order-1 (must have explicit row when m_order > 1)
            // - Second group: variables at time m_order (row can be implicit or explicit)

            // This pattern will match 
            std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

            auto iter = std::sregex_iterator(motif.begin(), motif.end(), pattern);

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

                        if (current_location >= bracket_ranges[1].first)
                            y_row = m_order;
                        else
                            throw std::logic_error("LHS of transition must specify time when m_order > 1");

                    } else
                        y_row = std::stoul(tmp_str);

                } else {

                    // If missing, we replace with the location
                    if (tmp_str != "")
                        y_row = std::stoul(tmp_str);
                    else
                        y_row = (current_location < bracket_ranges[1].first ? 0u: 1u);

                }

                if (y_row > m_order)
                    throw std::logic_error("The proposed row is out of range.");

                // Only variables at time m_order can be in the RHS (second bracketed group)
                if ((current_location >= bracket_ranges[1].first) && (y_row != m_order))
                    throw std::logic_error(
                        "Only the row " + std::to_string(m_order) +
                        " can be specified at the RHS of the motif."
                        );

                // And only variables before time m_order can be in the LHS
                if ((current_location < bracket_ranges[1].first) && (y_row >= m_order))
                    throw std::logic_error(
                        "Only rows 0 to " + std::to_string(m_order - 1) +
                        " can be specified at the LHS of the motif."
                        );

                if (selected[y_col * (m_order + 1) + y_row])
                    throw std::logic_error(
                        "The term " + i->str() + " shows more than once in the formula.");

                selected[y_col * (m_order + 1) + y_row] = true;

                locations.push_back(y_col * (m_order + 1) + y_row);
                signs.push_back(is_positive);
                

            }
        }
        else if (num_groups == m_order + 1)
        {
            // New behavior: each group corresponds to a time point (0 to m_order)
            
            // This pattern will match 
            std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

            auto iter = std::sregex_iterator(motif.begin(), motif.end(), pattern);

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

                // Determine which bracketed group this variable belongs to
                // Note: The regex pattern ensures all variables are within bracketed groups,
                // as the pattern only matches content inside brackets. Variables in covariate
                // expressions (after 'x') are handled separately by pattern_conditional above.
                size_t group_idx = 0;
                bool found_group = false;
                for (size_t g = 0; g < bracket_ranges.size(); ++g)
                {
                    if (current_location >= bracket_ranges[g].first && 
                        current_location < bracket_ranges[g].second)
                    {
                        group_idx = g;
                        found_group = true;
                        break;
                    }
                }
                
                // Safety check: ensure the variable was found in a bracketed group
                // This should never happen given the regex, but verify for safety
                if (!found_group)
                    throw std::logic_error(
                        "Internal error: variable " + i->str() + 
                        " not found within any bracketed group.");

                // Time location
                size_t y_row;
                std::string tmp_str = i->operator[](4u).str();
                
                // If row is explicitly specified, use it; otherwise infer from group position
                if (tmp_str != "")
                {
                    y_row = std::stoul(tmp_str);
                    
                    // Validate that explicit row matches the expected group position
                    if (y_row != group_idx)
                        throw std::logic_error(
                            "Explicit row index " + std::to_string(y_row) + 
                            " does not match the position in the formula (group " + 
                            std::to_string(group_idx) + ").");
                }
                else
                {
                    // Infer row from bracketed group position
                    y_row = group_idx;
                }

                if (y_row > m_order)
                    throw std::logic_error("The proposed row is out of range.");

                if (selected[y_col * (m_order + 1) + y_row])
                    throw std::logic_error(
                        "The term " + i->str() + " shows more than once in the formula.");

                selected[y_col * (m_order + 1) + y_row] = true;

                locations.push_back(y_col * (m_order + 1) + y_row);
                signs.push_back(is_positive);
                

            }
        }
        else
        {
            throw std::logic_error(
                "For a Markov model of order " + std::to_string(m_order) + 
                ", transition formulas must have either 2 bracketed groups " +
                "(for backwards compatibility) or exactly " + std::to_string(m_order + 1) + 
                " bracketed groups (found " + std::to_string(num_groups) + ").");
        }

        return;

    } 
    
    std::regex_match(formula, match, pattern_intercept);
    if (!match.empty())
    {

        extract_parts(match);

        // This pattern will match 
        std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

        auto iter = std::sregex_iterator(motif.begin(), motif.end(), pattern);

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

