import React, { useEffect, useState } from 'react';
import { queryTeamJobs } from 'shared/services/mock/jobs';
import JobList from 'components/jobs/JobList';

function TeamJobList({ teamId, forPart }) { // forPart can be "geneMapper" or "geneCruncher"
  const [jobs, setJobs] = useState([]);
  const [loading, setLoading] = useState(true);

  useEffect(async () => {
    setLoading(true);
    setJobs(await queryTeamJobs(teamId, forPart));
    setLoading(false);
  }, [teamId, forPart]);

  return <JobList isLoading={loading} jobs={jobs} />;
}

export default TeamJobList;
