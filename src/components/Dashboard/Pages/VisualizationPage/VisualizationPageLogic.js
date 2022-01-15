import { BACKEND_ADDRESS } from '../../../common/constants';
import { getAuthAndJsonHeader } from '../../../common/utils';

export default async function queryJob(id) {
  return fetch(`${BACKEND_ADDRESS}/job/${id}`, {
    headers: getAuthAndJsonHeader(),
  }).then((response) => {
    if (response.status !== 200) {
      throw Error(`Couldn't fetch job ${id}!`);
    }
    return response.json();
  });
}
