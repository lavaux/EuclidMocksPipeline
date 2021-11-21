from . import register_tool
from ..config import readConfig

@register_tool
def createProjectDirectoryTree(project: str = "Repo") -> None:
    """Setup a new repository of mocks and data files for a systematic study project

    Args:
        project (str, optional): [description]. Defaults to "Repo".
    """
    import os


    def project_mkdir(*subprojects):
        dir=os.path.join(project, *subprojects)
        print(f"  creating {dir}")
        os.mkdir(dir)


    print(f"Project: {project}")

    subdirs = [
        "Selections",
        "Catalogs4LE3",
        "GalaxyCatalogs",
        "Plots",
        "Cls",
        "NumberCounts",
        "RandomCatalogs",
        "2PCF",
        ("2PCF", "Measures"),
        ("2PCF", "Params"),
        ("2PCF", "Scripts"),
        "PK",
        ("PK", "Measures"),
        ("PK", "Params"),
        ("PK", "Scripts"),
        "CM-PK",
        ("CM-PK", "Measures"),
        ("CM-PK", "Params"),
        ("CM-PK", "Scripts"),
        "3PCF",
        ("3PCF", "Measures"),
        ("3PCF", "Params"),
        ("3PCF", "Scripts"),
        "BK",
        ("BK", "Measures"),
        ("BK", "Params"),
        ("BK", "Scripts")
    ]

    for s in subdirs:
        if type(s) == str:
            project_mkdir(s)
        else:
            project_mkdir(*s)

    # TODO: This is extremely bad! Fix that copy to something sane
    os.system(f"cp ../Repo/SelectionInputs/sel_input* ../{project}/Selections")