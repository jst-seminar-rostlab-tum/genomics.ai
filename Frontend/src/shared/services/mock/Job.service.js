const JobService = {
  async getJobs() {
    // TODO: add mock data
    return [];
  },

  async getTeamJobs(teamId, forPart) {
    const allJobs = await JobService.getJobs();
    // pretty random, this is just mocking
    if (forPart === 'geneMapper') {
      return allJobs.slice(teamId, allJobs.length / 2);
    }
    return allJobs.slice(allJobs.length / 2, allJobs.length - teamId);
  },
};

export default JobService;
