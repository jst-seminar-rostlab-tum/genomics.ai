import { BACKEND_ADDRESS } from '../../common/constants';
import { getAuthAndJsonHeader } from '../../common/utils';

async function queryJobs() {
  return fetch(`${BACKEND_ADDRESS}/jobs`, {
    headers: getAuthAndJsonHeader(),
  }).then((response) => {
    if (response.status !== 200) {
      throw Error("Couldn't fetch jobs");
    }
    return response.json();
  });
}
export default queryJobs;
