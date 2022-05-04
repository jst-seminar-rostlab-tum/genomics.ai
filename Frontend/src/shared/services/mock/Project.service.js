const mockProjects = [
  {
    _id: 'f7f7',
    forPart: 'geneMapper',
    status: 'DONE',
  },
  {
    _id: 'f3f6',
    forPart: 'geneCruncher',
    status: 'DONE',
  },
];

const ProjectService = {
  async getProjects() {
    return mockProjects;
  },

  async getTeamProjects(teamId, forPart) {
    const allProjects = await ProjectService.getProjects();
    return allProjects.filter((elem) => elem.forPart === forPart);
  },
};

export default ProjectService;
