const mockJobs = [
  {
    id: 'f7f7',
    forPart: 'geneMapper',
    status: 'DONE',
  },
  {
    id: 'f3f6',
    forPart: 'geneCruncher',
    status: 'DONE',
  },
];

const JobService = {
  async getJobs() {
    return mockJobs;
  },

  async getTeamJobs(teamId, forPart) {
    const allJobs = await JobService.getJobs();
    return allJobs.filter((elem) => elem.forPart === forPart);
  },
};

export default JobService;
