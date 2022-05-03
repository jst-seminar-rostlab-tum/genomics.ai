import React, { useEffect, useState } from 'react';
import JobService from 'shared/services/Job.service';
import JobList from 'components/jobs/JobList';

function TeamJobList({ team, forPart }) { // forPart can be "geneMapper" or "geneCruncher"
  const [jobs, setJobs] = useState([]);
  const [loading, setLoading] = useState(true);

  useEffect(async () => {
    setLoading(true);
    setJobs(await JobService.getTeamJobs(team.id, forPart));
    setLoading(false);
  }, [team, forPart]);

  return <JobList isLoading={loading} jobs={jobs} />;
}

export default TeamJobList;
