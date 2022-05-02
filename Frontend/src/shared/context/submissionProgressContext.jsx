import React, { useState } from 'react';
import { MULTIPART_UPLOAD_STATUS } from 'shared/utils/common/constants';

const SubmissionProgressContext = React.createContext();
SubmissionProgressContext.displayName = 'SubmissionProgressContext';

function SubmissionProgressProvider(props) {
  const [submissionProgress, setSubmissionProgress] = useState({
    status: MULTIPART_UPLOAD_STATUS.IDLE,
    uploadId: '',
    chunks: 0,
    uploaded: 0,
    remaining: [],
    uploadedParts: [],
  });

  const value = [submissionProgress, setSubmissionProgress];
  return <SubmissionProgressContext.Provider value={value} {...props} />;
}

function useSubmissionProgress() {
  const context = React.useContext(SubmissionProgressContext);
  if (context === undefined) {
    throw new Error('useSubmissionProgress must be used within a SubmissionProgressProvider');
  }
  return context;
}

export { SubmissionProgressProvider, useSubmissionProgress };