import React, { useState } from 'react';
import { MULTIPART_UPLOAD_STATUS } from 'shared/utils/common/constants';

const SubmissionProgressContext = React.createContext();
SubmissionProgressContext.displayName = 'SubmissionProgressContext';

const initSubmissionProgress = (uploadId) => ({
  status: MULTIPART_UPLOAD_STATUS.IDLE,
  uploadId,
  chunks: 0,
  uploaded: 0,
  remaining: [],
  uploadedParts: [],
});

function SubmissionProgressProvider(props) {
  const [submissionProgress, setSubmissionProgress] = useState({});

  const value = [submissionProgress, setSubmissionProgress];

  // eslint-disable-next-line react/jsx-props-no-spreading
  return <SubmissionProgressContext.Provider value={value} {...props} />;
}

function useSubmissionProgress() {
  const context = React.useContext(SubmissionProgressContext);
  if (context === undefined) {
    throw new Error('useSubmissionProgress must be used within a SubmissionProgressProvider');
  }
  return context;
}

export { SubmissionProgressProvider, useSubmissionProgress, initSubmissionProgress };
