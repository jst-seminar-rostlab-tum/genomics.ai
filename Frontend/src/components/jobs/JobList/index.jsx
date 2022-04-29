import React from 'react';
import LoadingList from 'components/general/LoadingList';
import JobCard from 'components/jobs/JobCard';

function JobList({ isLoading, jobs }) {
  return (
    <LoadingList
      isLoading={isLoading}
      elements={jobs}
      cardBuilder={(job) => <JobCard job={job} />}
      noElementsMessage="No jobs."
    />
  );
}

export default JobList;
