import { BACKEND_ADDRESS } from '../../utils/common/constants';
import { getAuthAndJsonHeader } from '../../utils/common/utils';

export default async function queryJobs() {
  return fetch(`${BACKEND_ADDRESS}/jobs`, {
    headers: getAuthAndJsonHeader(),
  }).then((response) => {
    if (response.status !== 200) {
      throw Error("Couldn't fetch jobs");
    }
    return response.json();
  }).then((result) => result.map((job) => ({ id: job._id, ...job })));
}

// forPart can be "geneMapper" or "geneCruncher"
export async function queryTeamJobs(teamId, forPart) {
  const allJobs = await queryJobs();
  // pretty random, this is just mocking
  if (forPart === 'geneMapper') {
    return allJobs.slice(teamId, allJobs.length / 2);
  }
  return allJobs.slice(allJobs.length / 2, allJobs.length - teamId);
}
