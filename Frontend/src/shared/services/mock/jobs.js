import { BACKEND_ADDRESS } from '../../utils/common/constants';
import { getAuthAndJsonHeader } from '../../utils/common/utils';

let mockJobs = [
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
  return mockJobs.filter((elem) => elem.forPart === forPart);
}
