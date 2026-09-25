import cobra


def two_route_model():
    """A -> B by one step that draws expensive enzyme, or by two steps that draw none."""
    model = cobra.Model("routes")
    a, b, c, enzyme = (cobra.Metabolite(x) for x in "ABCE")

    def reaction(name, stoichiometry, upper=1000.0):
        r = cobra.Reaction(name, lower_bound=0.0, upper_bound=upper)
        r.add_metabolites(stoichiometry)
        model.add_reactions([r])

    reaction("supply", {a: 1}, upper=1.0)
    reaction("direct", {a: -1, b: 1, enzyme: -10})
    reaction("usage_prot_E", {enzyme: 1})
    reaction("step1", {a: -1, c: 1})
    reaction("step2", {c: -1, b: 1})
    reaction("demand", {b: -1})
    model.objective = "demand"
    return model

