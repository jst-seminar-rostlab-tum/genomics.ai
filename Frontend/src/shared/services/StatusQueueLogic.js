import { BACKEND_ADDRESS } from '../utils/common/constants';
import { getAuthAndJsonHeader } from '../utils/common/utils';

export async function queryJobs() {
  return fetch(`${BACKEND_ADDRESS}/jobs`, {
    headers: getAuthAndJsonHeader(),
  }).then((response) => {
    if (response.status !== 200) {
      throw Error("Couldn't fetch jobs");
    }
    return response.json();
  });
}

export function filterInProgress(jobs) {
  return jobs.filter((job) => job.status === 'DONE');
}
