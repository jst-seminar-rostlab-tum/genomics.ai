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

  getOwnProjects: async () => ProjectService.getProjects(),

  async getProject(id) {
    const project = mockProjects.find((p) => String(p._id) === String(id));
    return project ?? {
      _id: id,
      owner: 'mock',
      uploadId: '12345-abcde',
      atlasId: '626ea3311d7d1a27de465b63',
      fileName: 'test1.h5ad',
      location: `./testData/test_file${id}.csv`,
      fileSize: 10000000,
      modelId: '626ea2361d7d1a27de465b5e',
      name: `Test project ${id}`,
      resultSize: '10000',
      status: 'DONE',
      uploadDate: '2022-04-30T106:00:00.000Z',
    };
  },
};

export default ProjectService;
