// by https://www.npmjs.com/package/react-image-crop
import { useEffect } from 'react';

export default function useDebounceEffect(
  fn,
  waitTime,
  deps,
) {
  useEffect(() => {
    const t = setTimeout(() => {
      // eslint-disable-next-line prefer-spread
      fn.apply(undefined, deps);
    }, waitTime);

    return () => {
      clearTimeout(t);
    };
  }, deps);
}
