import React, { useEffect, useState } from 'react';
import { queryTeamJobs } from 'shared/services/mock/jobs';
import JobList from 'components/jobs/JobList';

function TeamJobList({ team, forPart }) { // forPart can be "geneMapper" or "geneCruncher"
  const [jobs, setJobs] = useState([]);
  const [loading, setLoading] = useState(true);

  useEffect(async () => {
    setLoading(true);
    setJobs(await queryTeamJobs(team.id, forPart));
    setLoading(false);
  }, [team, forPart]);

  return <JobList isLoading={loading} jobs={jobs} />;
}

export default TeamJobList;
